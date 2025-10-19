#!/usr/bin/env python3
import sys
import yaml
import demes
import subprocess
import argparse
import re
import random
from collections import defaultdict, deque

# ---------- File paths ----------
OTHER_YAML      = "src/other.yaml"
DEMES_YAML      = "src/demes.yaml"
EXECUTABLE      = "./executables/labp_v19"

# ---------- Defaults ----------
parameters = {
    "seed":         "1",
    "nruns":        "1",
    "kingman_coal": "1",
    "drift_sim":    "0",
    "msOutput":     "1",
    "popSizeVec":   "1000 1000",
    "inv_freq":     "0.2 0.3",
    "originPop":    "pop1",                  # 0-based index or name of origin deme (e.g. "pop2")
    "inv_size":     "1000000",            # inversion length (bp)
    "speciation":   "0 0 0",
    "demography":   "0 0 0",
    "inv_age":      "-1",                 # −1 = no inversion; ≥0 = generation where inversion originates
    "migRate":      "0.02",
    "BasesPerMorgan":"1e8",
    "randPhi":      "0",
    "phi":          "0.2",
    "invRange":     "0 1e3",
    "fixedSNPs":    "1 10 1e-6",
    "n_SNPs":       "0",
    "snpPositions": "500 550 600 650 700 750 800 850 900 950",
    "randomSample": "0",
    "tempRead":     "10 0",
}

# ---------- Helpers for demes ----------
def epoch_start_size(epoch):
    return int(epoch.start_size)

def leaf_demes(graph):
    return [d for d in graph.demes if d.epochs and d.epochs[-1].end_time == 0]

def pop_index_order(leaves):
    def keyfn(name):
        m = re.fullmatch(r"pop(\d+)", name)
        if m: return (0, int(m.group(1)))
        return (1, name)
    names = sorted([d.name for d in leaves], key=keyfn)
    return {name: i for i, name in enumerate(names)}

def build_children_map(graph):
    from collections import defaultdict
    children = defaultdict(list)
    for d in graph.demes:
        if getattr(d, "ancestors", None):
            if len(d.ancestors) > 1:
                raise NotImplementedError(
                    f"Admixture is not supported for this parser (deme={d.name}, ancestors={d.ancestors})."
                )
            parent = d.ancestors[0]
            children[parent].append(d.name)
    for p in children:
        children[p].sort(key=lambda x: (x.startswith("pop"), x))
    return children

def leaves_under(node, children, leafset):
    out = []
    from collections import deque
    q = deque([node])
    while q:
        u = q.popleft()
        if u in leafset:
            out.append(u)
        for v in children.get(u, []):
            q.append(v)
    return out

def find_roots(graph):
    return [d.name for d in graph.demes if len(getattr(d, "ancestors", [])) == 0]

def pick_root(graph):
    roots = find_roots(graph)
    if not roots:
        raise ValueError("No root deme (a deme with no ancestors) found in demes.yaml")
    if len(roots) == 1:
        return roots[0]
    by_name = {d.name: d for d in graph.demes}
    roots.sort(key=lambda n: by_name[n].start_time, reverse=True)
    return roots[0]

def fmt_time(x):
    try:
        if x == float("inf"):
            return "inf"
        return f"{x:g}"
    except Exception:
        return str(x)

def ancestor_size_at(deme, T):
    eps = 1e-9
    for ep in deme.epochs:
        if abs(ep.end_time - T) < eps:
            return int(ep.start_size)
    for ep in deme.epochs:
        if ep.end_time > T:
            return int(ep.start_size)
    return int(deme.epochs[-1].start_size)

def render_tree_ascii(graph, only_these_leaves=None):
    children = build_children_map(graph)
    root = pick_root(graph)
    by_name = {d.name: d for d in graph.demes}

    if only_these_leaves:
        needed = set()
        parents = {c: (by_name[c].ancestors[0] if getattr(by_name[c], "ancestors", None) else None)
                   for c in by_name}
        for leaf in only_these_leaves:
            u = leaf
            while u is not None and u not in needed:
                needed.add(u)
                u = parents[u]
        pruned = {p: [c for c in cs if c in needed] for p, cs in children.items() if p in needed}
        children = pruned

    lines = []
    def dfs(u, prefix="", depth=0, is_last=True):
        connector = "" if depth == 0 else ("└─ " if is_last else "├─ ")
        label = f"{u}  (start={fmt_time(by_name[u].start_time)})"
        lines.append(prefix + connector + label)
        kids = children.get(u, [])
        for i, v in enumerate(kids):
            last = (i == len(kids) - 1)
            next_prefix = prefix if depth == 0 else prefix + ("   " if is_last else "│  ")
            dfs(v, next_prefix, depth + 1, last)
    dfs(root)
    return "\n".join(lines)

# ---------- Speciation (sizes + events) ----------
def build_speciation_from_demes(graph, ancestor_freqs=None, default_F=0.2):
    leaves = leaf_demes(graph)
    if not leaves:
        raise ValueError("No present-day demes (end_time==0) found.")
    idx_map = pop_index_order(leaves)
    leaf_names_in_order = [name for name, _ in sorted(idx_map.items(), key=lambda kv: kv[1])]

    sizes = []
    for name in leaf_names_in_order:
        d = next(dd for dd in leaves if dd.name == name)
        sizes.append(str(epoch_start_size(d.epochs[-1])))

    children = build_children_map(graph)
    by_name = {d.name: d for d in graph.demes}
    leafset = set(leaf_names_in_order)
    active_groups = [[name] for name in leaf_names_in_order]

    def find_group_idx(leaf):
        for i, g in enumerate(active_groups):
            if leaf in g:
                return i
        raise RuntimeError(f"Leaf {leaf} not found in current groups")

    internal_nodes = list(children.keys())
    internal_nodes.sort(key=lambda n: min((by_name[c].start_time for c in children[n])), reverse=False)

    events = []
    plan_log = []
    groups_log = []

    def get_parent_F(parent):
        if ancestor_freqs is None:
            return float(default_F)
        if parent not in ancestor_freqs:
            raise ValueError(f"ancestor_frequencies is missing a value for internal node '{parent}'")
        F = float(ancestor_freqs[parent])
        if not (0.0 <= F <= 1.0):
            raise ValueError(f"ancestor_frequencies['{parent}'] must be in [0,1], got {F}")
        return F

    def leaves_under_name(name):
        return sorted(leaves_under(name, children, leafset))

    for parent in internal_nodes:
        child_sets = []
        for child in children[parent]:
            under = leaves_under_name(child)
            if under:
                child_sets.append((child, under))
        if len(child_sets) < 2:
            continue

        T = min(by_name[ch].start_time for ch, _ in child_sets)
        R = ancestor_size_at(by_name[parent], T)
        F_parent = get_parent_F(parent)

        child_group_idxs = []
        for ch, under in child_sets:
            gi = find_group_idx(under[0])
            child_group_idxs.append((gi, ch, under))

        child_group_idxs.sort(key=lambda x: x[0])
        sink_idx = child_group_idxs[0][0]
        sink_child = child_group_idxs[0][1]

        for gi, ch, under in child_group_idxs[1:]:
            A = sink_idx
            B = gi
            events.append((A, B, T, F_parent, R))
            plan_log.append(
                f"t={fmt_time(T)}: merge group {B} ({ch}:{under}) → group {A} ({sink_child}); F={F_parent:g}, R={R}"
            )
            active_groups[A].extend(active_groups[B])
            del active_groups[B]
            groups_log.append([list(g) for g in active_groups])

    pop_sizes_str = " ".join(sizes)
    if events:
        speciation_flat = ["1"]
        for (A, B, T, F, R) in events:
            speciation_flat += [str(A), str(B), f"{T:g}", f"{F:g}", str(R)]
        speciation_str = " ".join(speciation_flat)
    else:
        speciation_str = "0 0 0"

    debug_info = {
        "leaf_order": leaf_names_in_order,
        "initial_groups": [[x] for x in leaf_names_in_order],
        "events": events,
        "plan_log": plan_log,
        "groups_log": groups_log,
    }
    return pop_sizes_str, speciation_str, debug_info

# ---------- Demography (built from all demes) ----------
def build_demography_from_demes_full(graph, leaf_order, spec_events):
    by_name = {d.name: d for d in graph.demes}
    npops   = len(leaf_order)

    children = build_children_map(graph)
    leafset  = set(leaf_order)

    def leaves_under_name(name):
        return sorted(leaves_under(name, children, leafset))

    demog_records = defaultdict(list)
    for d in graph.demes:
        for ep in d.epochs:
            if ep.end_time and ep.end_time > 0:
                demog_records[ep.end_time].append((d.name, int(ep.start_size)))

    if not demog_records:
        return "0 0 0", []

    merges_at_time = defaultdict(list)
    for (A, B, T, _, _) in spec_events:
        merges_at_time[T].append((A, B))

    active_groups = [[name] for name in leaf_order]

    def group_index_for_deme(name):
        target = set(leaves_under_name(name))
        for i, g in enumerate(active_groups):
            if set(g) == target:
                return i
        return None

    demog_entries = ["1"]
    dbg_lines = []

    for t in sorted(demog_records.keys()):
        sizes = [0] * npops
        applied = []

        for name, sz in demog_records[t]:
            gi = group_index_for_deme(name)
            if gi is not None and gi < len(active_groups):
                sizes[gi] = sz
                applied.append((gi, name, sz))

        if any(v != 0 for v in sizes[:len(active_groups)]):
            demog_entries.append(f"{t:g}")
            demog_entries += [str(v) for v in sizes]
            applied_str = ", ".join([f"g{gi}({name})={sz}" for gi, name, sz in applied])
            dbg_lines.append(f"  t={fmt_time(t)}: {applied_str}")

        for (A, B) in sorted(merges_at_time.get(t, []), key=lambda x: (x[0], x[1])):
            if A < len(active_groups) and B < len(active_groups):
                active_groups[A].extend(active_groups[B])
                del active_groups[B]

    return " ".join(demog_entries), dbg_lines

# ---------- Load other.yaml ----------
try:
    with open(OTHER_YAML, 'r') as f:
        other_params = yaml.safe_load(f) or {}
except Exception as e:
    print(f"Error loading {OTHER_YAML}: {e}", file=sys.stderr)
    sys.exit(1)

# Copy keys from other.yaml
for key, val in (other_params or {}).items():
    # copy any known key from defaults
    if key in parameters:
        parameters[key] = str(val)
        continue
    # ALSO copy any per-pop sample strings like nCarriers, nCarriers1, ...
    if key.startswith("nCarriers"):
        parameters[key] = str(val)

ancestor_freqs = other_params.get("ancestor_frequencies", None)
if ancestor_freqs is not None and not isinstance(ancestor_freqs, dict):
    print("Error: 'ancestor_frequencies' in other.yaml must be a mapping (dict).", file=sys.stderr)
    sys.exit(1)

DEFAULT_ANCESTOR_F = 0.2

# ---------- Load demes.yaml & derive sizes/speciation/migration/demography ----------
try:
    graph = demes.load(DEMES_YAML)
except Exception as e:
    print(f"Error loading {DEMES_YAML}: {e}", file=sys.stderr)
    sys.exit(1)

try:
    pop_sizes_str, speciation_str, debug_info = build_speciation_from_demes(
        graph,
        ancestor_freqs=ancestor_freqs,
        default_F=DEFAULT_ANCESTOR_F,
    )
    parameters["popSizeVec"] = pop_sizes_str
    parameters["speciation"] = speciation_str
except (NotImplementedError, ValueError) as e:
    print(f"[speciation parser] {e}", file=sys.stderr)
    sys.exit(1)

demog_str, demog_dbg = build_demography_from_demes_full(
    graph,
    debug_info["leaf_order"],
    debug_info["events"]
)
parameters["demography"] = demog_str

# --------------------------------------------------------------------
# Inversion configuration (present-day inv_freq + originPop at inv_age)
#   - inv_freq: present-day per-deme frequencies (0..N non-zero allowed)
#   - originPop: REQUIRED iff inv_age >= 0 (index "2" or leaf name "pop2")
#   - inv_age == -1 disables the boundary (ignore originPop)
#   - This uses present-day leaf order from debug_info["leaf_order"] (C++ uses same order)
#   - Computes origin deme's start time from demes.yaml (C++ will re-infer it)
# --------------------------------------------------------------------
try:
    inv_age_val = int(parameters.get("inv_age", "-1") or "-1")
except ValueError:
    inv_age_val = -1  # -1 = no inversion-origin event

leaf_names = debug_info["leaf_order"]  # present-day demes in C++ order
pops = len(leaf_names)

# Validate inv_freq values
raw_inv_freq = (parameters.get("inv_freq", "") or "").strip()
inv_freqs = [float(x) for x in raw_inv_freq.split()] if raw_inv_freq else []
if inv_freqs:
    if len(inv_freqs) != pops:
        print(f"[Error] inv_freq length ({len(inv_freqs)}) != number of present-day demes ({pops}).",
              file=sys.stderr)
        sys.exit(1)
    bad = [(i, f) for i, f in enumerate(inv_freqs) if not (0.0 <= f <= 1.0)]
    if bad:
        bad_str = ", ".join([f"leaf {i}({leaf_names[i]}): {f}" for i, f in bad])
        print(f"[Error] inv_freq entries must be in [0,1]. Offenders: {bad_str}", file=sys.stderr)
        sys.exit(1)

if inv_age_val < 0:
    # Feature OFF → ignore originPop if user provided it
    if "originPop" in parameters or "originPop" in (other_params or {}):
        print("[Warning] inv_age == -1; ignoring originPop (no inversion-origin collapse).",
              file=sys.stderr)
else:
    # inv_age >= 0 → originPop is REQUIRED (index or leaf name)
    origin_value = (other_params or {}).get("originPop", "").strip() or parameters.get("originPop", "").strip()
    if not origin_value:
        print("[Error] inv_age >= 0 requires 'originPop' (index or present-day leaf name) in other.yaml.",
              file=sys.stderr)
        sys.exit(1)

    # Accept "2" or "pop2" and map to the present-day leaf index
    if origin_value.isdigit():
        origin_idx = int(origin_value)
        if not (0 <= origin_idx < pops):
            print(f"[Error] originPop index {origin_idx} out of range for {pops} leaves.",
                  file=sys.stderr)
            sys.exit(1)
        origin_name = leaf_names[origin_idx]
    else:
        try:
            origin_idx = leaf_names.index(origin_value)
            origin_name = origin_value
        except ValueError:
            print(f"[Error] originPop name '{origin_value}' not found among present-day leaves {leaf_names}.",
                  file=sys.stderr)
            sys.exit(1)

    parameters["originPop"] = str(origin_idx)
    print(f"[Info] originPop set to {origin_idx} ({origin_name}). "
          f"At inv_age={inv_age_val}, all inverted lineages collapse here and convert to standard.")

    # --------------------------------------------------------------------
    # Compute origin_start_time (your definition: use the deme's end_time)
    # --------------------------------------------------------------------
    if inv_age_val >= 0:
        origin_name = leaf_names[origin_idx]  # present-day leaf name for originPop
        by_name = {d.name: d for d in graph.demes}
        if origin_name not in by_name:
            print(f"[Error] originPop '{origin_name}' not found in demes.yaml.", file=sys.stderr)
            sys.exit(1)

        origin_deme = by_name[origin_name]
        origin_start_time = float(origin_deme.end_time)

        # Optional: sanity checks
        if origin_start_time > inv_age_val:
            print(f"[Error] origin_start_time ({origin_start_time:g}) > inv_age ({inv_age_val:g}).",
                file=sys.stderr)
            sys.exit(1)

        print(f"[Info] origin_start_time={origin_start_time:g} "
            f"(taken from deme '{origin_name}' end_time in demes.yaml).")

# --------------------------------------------------------------------
# invRange handling (ALWAYS pass bp to C++; C++ converts using BasesPerMorgan)
# Precedence: explicit invRange in other.yaml > computed from inv_size
# --------------------------------------------------------------------
def _parse_two_floats(s: str):
    parts = [p for p in (s or "").split() if p.strip()]
    if len(parts) != 2:
        raise ValueError(f"invRange must have exactly two numbers, got: '{s}'")
    return float(parts[0]), float(parts[1])

explicit_invRange = (other_params or {}).get("invRange", None)
if explicit_invRange:
    try:
        L_bp, R_bp = _parse_two_floats(str(explicit_invRange))
        if L_bp > R_bp:
            print("[Error] invRange must be 'start end' with start <= end (in bp).", file=sys.stderr)
            sys.exit(1)
        parameters["invRange"] = f"{L_bp:g} {R_bp:g}"
        print(f"[Info] Using invRange from other.yaml (bp): {parameters['invRange']}")
    except Exception as e:
        print(f"[Error] Bad invRange in other.yaml: {e}", file=sys.stderr)
        sys.exit(1)
else:
    size_bp = float(parameters.get("inv_size", "0") or "0")
    if size_bp > 0:
        parameters["invRange"] = f"0 {size_bp:g}"
        print(f"[Info] invRange derived from inv_size (bp): {parameters['invRange']}")
    else:
        if inv_age_val >= 0:
            print("[Warning] inv_age ≥ 0 but neither invRange nor inv_size provided; "
                  "C++ will use default. Consider setting inv_size or invRange.",
                  file=sys.stderr)

def _extract_default_mig_rate_from_yaml(path):
    try:
        with open(path, "r") as _f:
            raw = yaml.safe_load(_f) or {}
        dflt = (raw.get("defaults") or {}).get("migration") or {}
        rate = dflt.get("rate", None)
        return float(rate) if rate is not None else None
    except Exception:
        return None

mig_rate = None
try:
    if graph.migrations and hasattr(graph.migrations[0], "rate"):
        mig_rate = float(graph.migrations[0].rate)
except Exception:
    mig_rate = None
if mig_rate is None:
    mig_rate = _extract_default_mig_rate_from_yaml(DEMES_YAML)
if mig_rate is not None:
    parameters["migRate"] = f"{mig_rate:g}"

# ---------- Debug prints ----------
print("\n=== Parsed speciation tree (past → present) ===")
ascii_tree = render_tree_ascii(graph, only_these_leaves=debug_info["leaf_order"])
print(ascii_tree)

print("\n=== Present-day leaves (order → indices) ===")
for i, name in enumerate(debug_info["leaf_order"]):
    print(f"  {i}: {name}")

print("\n=== Speciation plan (after reindexing between events) ===")
if debug_info["events"]:
    for line in debug_info["plan_log"]:
        print("  " + line)
    if debug_info["groups_log"]:
        print("\n--- Group reindex snapshots ---")
        for step, groups in enumerate(debug_info["groups_log"], 1):
            print(f"  after merge {step}: {groups}")
else:
    print("  (none)")

print("\n=== Demography plan (from demes.yaml) ===")
if parameters["demography"].startswith("1"):
    if demog_dbg:
        for line in demog_dbg:
            print(line)
    else:
        print("  (events exist, but all mapped to times where affected groups didn't yet exist pre-merge)")
else:
    print("  (none)")

# ---------- Allow user tweaks EXCEPT derived fields ----------
def print_parameters(params):
    for k, v in params.items():
        print(f"{k}: {v}")

def modify_params(parameters):
    print("Current Parameters:")
    print_parameters(parameters)
    ans = input("\nDo you want to modify any parameters? (Y/N): ").strip().lower()
    if ans == 'y':
        print("(Note: 'popSizeVec', 'speciation', and 'demography' are derived from demes.yaml and will be ignored if changed.)")
        user_input = input("\nEnter parameters to modify (e.g., --seed 2 --nruns 3): ")
        parser = argparse.ArgumentParser()
        for key in parameters:
            if key in ["popSizeVec", "speciation", "demography"]:
                continue
            if key in ["inv_freq","invRange","fixedSNPs","tempRead","nCarriers","snpPositions",
                       "nCarriers1","nCarriers2","nCarriers3","nCarriers4","nCarriers5","nCarriers6","nCarriers7","nCarriers8","nCarriers9"]:
                parser.add_argument(f"--{key}", nargs='+')
            else:
                parser.add_argument(f"--{key}")
        if user_input.strip():
            args = parser.parse_args(user_input.split())
            for key, val in vars(args).items():
                if val is not None:
                    parameters[key] = " ".join(val) if isinstance(val, list) else str(val)
        print("\nUpdated Parameters:")
        print_parameters(parameters)
    else:
        print("\nNo changes made.")

modify_params(parameters)

# Re-assert derived values (in case user tried to change them)
parameters["popSizeVec"] = pop_sizes_str
parameters["speciation"] = speciation_str
parameters["demography"] = demog_str

# ---------- Validation helpers for sample arguments ----------
def _count_pops_from_graph(graph) -> int:
    return len(leaf_demes(graph))

def _ints_in_str(s: str):
    s = (s or "").strip()
    if not s:
        return []
    try:
        return [int(x) for x in s.split()]
    except ValueError:
        raise ValueError(f"Expected integers in '{s}'.")

def _collect_per_pop_strings(parameters: dict, pops: int, random_flag: str):
    """
    Returns a list of sample strings to send at the tail of argv,
    matching the C++ expectations:
      - If random == "1": returns [ tempRead ] (must contain `pops` ints)
      - If random == "0": returns [ tempRead, nCarriers, nCarriers1, ... ] (length == pops),
                          each must contain exactly 2 ints.
    Raises ValueError on mismatch.
    """
    if random_flag == "1":
        if "tempRead" not in parameters:
            raise ValueError("randomSample==1 requires 'tempRead'.")
        vals = _ints_in_str(parameters["tempRead"])
        if len(vals) != pops:
            raise ValueError(f"randomSample==1: 'tempRead' must have {pops} integers, got {len(vals)}.")
        return [parameters["tempRead"]]

    # random_flag == "0"
    per_pop = []
    # pop0 comes from tempRead
    vals0 = _ints_in_str(parameters.get("tempRead", ""))
    if len(vals0) != 2:
        raise ValueError(f"randomSample==0: pop0 string 'tempRead' must have exactly 2 integers, got {len(vals0)}.")
    per_pop.append(parameters["tempRead"])

    # remaining pops from nCarriers, nCarriers1, nCarriers2, ...
    carriers_keys = ["nCarriers"] + [f"nCarriers{i}" for i in range(1, pops-1)]
    for idx, key in enumerate(carriers_keys, start=1):
        if key not in parameters:
            raise ValueError(f"randomSample==0: missing per-pop string for pop{idx}: '{key}'.")
        v = parameters[key]
        ints = _ints_in_str(v)
        if len(ints) != 2:
            raise ValueError(f"randomSample==0: '{key}' for pop{idx} must have exactly 2 integers, got {len(ints)}.")
        per_pop.append(v)

    if len(per_pop) != pops:
        raise ValueError(f"randomSample==0: expected {pops} per-pop strings, got {len(per_pop)}.")
    return per_pop

# ---------- Final assembly ----------
# Base (fixed-order) arguments up to and including randomSample:
base_keys_order = [
    "seed","nruns","kingman_coal","drift_sim","msOutput",
    "popSizeVec","inv_freq","originPop","inv_size","speciation","demography","inv_age",
    "migRate","BasesPerMorgan","randPhi","phi","invRange","fixedSNPs",
    "n_SNPs","snpPositions","randomSample",
]

# Build the tail according to randomSample rule
pops = _count_pops_from_graph(graph)
random_flag = parameters.get("randomSample","0").strip()
per_pop_strings = _collect_per_pop_strings(parameters, pops, random_flag)

# Compose final argv list
args_list = [parameters[k] for k in base_keys_order] + per_pop_strings

print("\nFinal parameters being passed:")
for k in base_keys_order:
    print(f"{k}: {parameters[k]}")
# Print the tail clearly
if random_flag == "1":
    print(f"tempRead: {parameters['tempRead']}")
else:
    print(f"tempRead (pop0): {parameters['tempRead']}")
    for i in range(1, pops):
        key = "nCarriers" if i == 1 else f"nCarriers{i-1}"
        print(f"{key} (pop{i}): {parameters[key]}")

try:
    result = subprocess.run([EXECUTABLE] + args_list, capture_output=True, text=True)
except FileNotFoundError:
    print(f"\nERROR: Executable not found at {EXECUTABLE}", file=sys.stderr)
    sys.exit(1)

print("\nC++ Program Output:\n", result.stdout)
if result.stderr:
    print("\nErrors:\n", result.stderr)
    with open("Output_log.txt", "w") as log_file:
        log_file.write(result.stderr)
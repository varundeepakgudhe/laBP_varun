#!/usr/bin/env python3
import sys
import yaml
import subprocess
import argparse
import json

# Default file names for the YAML files.
OTHER_YAML      = "src/other.yaml"
DEMOGRAPHY_JSON = "src/demo_dict.json"   # your sparse‐demography map

# The default parameters dictionary
parameters = {
    "seed":         "1",
    "nruns":        "1",
    "kingman_coal": "1",
    "drift_sim":    "0",
    "msOutput":     "1",
    "popSizeVec":   "10000 10000 10000",
    "inv_freq":     "0.2 0.4 0.1",
    "speciation":   "1 0 2 1000 0.3",
    "demography":   "0 0 0",
    "inv_age":      "0",
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
    "nCarriers":    "10 0"
}

# 1) Load other.yaml and override defaults
try:
    with open(OTHER_YAML, 'r') as f:
        other_params = yaml.safe_load(f)
except Exception as e:
    print(f"Error loading {OTHER_YAML}: {e}", file=sys.stderr)
    sys.exit(1)

# Override *all* keys from other.yaml
for key, val in other_params.items():
    if key in parameters:
        parameters[key] = str(val)

# 2) (Optional) Interactive tweak
def print_parameters(params):
    for k, v in params.items():
        print(f"{k}: {v}")

def modify_params(parameters):
    print("Current Parameters:")
    print_parameters(parameters)
    if input("\nDo you want to modify any parameters? (Y/N): ").lower() == 'y':
        user_input = input("\nEnter parameters to modify (e.g., --seed 2 --nruns 3): ")
        parser = argparse.ArgumentParser()
        for key in parameters:
            # list‐valued vs scalar
            parser.add_argument(f"--{key}", nargs='+' if ' ' in parameters[key] else None)
        args = parser.parse_args(user_input.split())
        for key, val in vars(args).items():
            if val is not None:
                parameters[key] = " ".join(val) if isinstance(val, list) else str(val)
        print("\nUpdated Parameters:")
        print_parameters(parameters)
    else:
        print("\nNo changes made.")

modify_params(parameters)

# 3) (Optional) Load your demography JSON map for explicit demography events
try:
    with open(DEMOGRAPHY_JSON, 'r') as f:
        demog_map = json.load(f)
except FileNotFoundError:
    demog_map = {}

if demog_map:
    demog_entries = ["1"]
    for t_str, sizes in sorted(demog_map.items(), key=lambda kv: float(kv[0])):
        demog_entries.append(t_str)
        demog_entries += [str(int(x)) for x in sizes]
    parameters["demography"] = " ".join(demog_entries)

# 4) Emit final list in the exact order main.cpp expects
keys_order = [
    "seed","nruns","kingman_coal","drift_sim","msOutput",
    "popSizeVec","inv_freq","speciation","demography","inv_age",
    "migRate","BasesPerMorgan","randPhi","phi","invRange","fixedSNPs",
    "n_SNPs","snpPositions","randomSample","tempRead","nCarriers"
]
args_list = [parameters[k] for k in keys_order]

print("\nFinal parameters being passed:")
for k in keys_order:
    print(f"{k}: {parameters[k]}")

# 5) Invoke the C++ executable
command = ["./executables/labp_v19"] + args_list
result = subprocess.run(command, capture_output=True, text=True)

print("\nC++ Program Output:\n", result.stdout)
if result.stderr:
    print("\nErrors:\n", result.stderr)
    with open("Output_log.txt", "w") as log_file:
        log_file.write(result.stderr)

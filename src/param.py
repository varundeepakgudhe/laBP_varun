#!/usr/bin/env python3
import sys
import yaml
import demes
import subprocess

parameters = {
    "seed": "-1",
    "nruns": "1",
    "kingman_coal": "1",
    "drift_sim": "0",
    "msOutput": "1",
    "popSizeVec": "10000 10000",
    "inv_freq": "0 0",
    "speciation": "1 10000 0",
    "demography": "0 0 0",
    "inv_age": "0",
    "migRate": "0",
    "BasesPerMorgan": "1e8",
    "randPhi": "0",
    "phi": "0",
    "invRange": "0 1e3",
    "fixedSNPs": "1 10 1e-6",
    "n_SNPs": "1",
    "snpPositions": "998 999",
    "randomSample": "0",
    "tempRead": "10 0",
    "nCarriers": "10 0"
}

if len(sys.argv) < 3:
    print("Usage: python param.py <other.yaml> <demes.yaml>")
    sys.exit(1)

with open(sys.argv[1], 'r') as f:
    other_params = yaml.safe_load(f)

for key in other_params:
    if key in parameters and key not in ["popSizeVec", "speciation", "migRate"]:
        parameters[key] = str(other_params[key])

graph = demes.load(sys.argv[2])

pop_sizes = []
for deme in graph.demes:
    if deme.name != "ancestor":
        pop_sizes.append(str(deme.epochs[0].start_size))
parameters["popSizeVec"] = " ".join(pop_sizes)

ancestor_end_time = None
for deme in graph.demes:
    if deme.name == "ancestor":
        ancestor_end_time = deme.epochs[0].end_time
        break
if ancestor_end_time is not None:
    speciation_parts = parameters["speciation"].split()
    if len(speciation_parts) >= 3:
        speciation_parts[1] = str(ancestor_end_time)
        parameters["speciation"] = " ".join(speciation_parts)

if graph.migrations:
    parameters["migRate"] = str(graph.migrations[0].rate)

keys_order = ["seed", "nruns", "kingman_coal", "drift_sim", "msOutput",
              "popSizeVec", "inv_freq", "speciation", "demography", "inv_age",
              "migRate", "BasesPerMorgan", "randPhi", "phi", "invRange", "fixedSNPs",
              "n_SNPs", "snpPositions", "randomSample", "tempRead", "nCarriers"]

args = [parameters[key] for key in keys_order]
print("Updated parameters:")
for key in keys_order:
    print(f"{key}: {parameters[key]}")

command = ["./executables/labp_v19"] + args
result = subprocess.run(command, capture_output=True, text=True)

print("C++ Program Output:\n", result.stdout)
if result.stderr:
    print("Errors:\n", result.stderr)
    with open("Output_log.txt", "w") as log_file:
        log_file.write(result.stderr)
import subprocess
import argparse

def run_cpp_program(params):
    # Call the compiled C++ program and pass the parameters
    command = ["./executables/labp_v19"] + params
    result = subprocess.run( command, capture_output=True, text=True)

    # Print the output received from the C++ program
    print("C++ Program Output:\n", result.stdout)
    if result.stderr:
        print("Errors:\n", result.stderr)
        with open("Output_log.txt", "w") as log_file:
            log_file.write(result.stderr)

## Default parameters
parameters = {
    "seed": "-1",
    "nruns": "10",
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
    "fixedSNPs": "0 1000 1e-6",
    "n_SNPs": "1",
    "snpPositions": "1000.01 1000.011",
    "randomSample": "0",
    "tempRead": "10 0",
    "nCarriers": "10 0"
}

## Display all the default parameters.
def print_parameters(parameters):
   
    for key, value in parameters.items():
        print(f"{key}: {value}")


## Modify/Update parameters.
def modify_params(parameters):
    print("Current Parameters:")
    print_parameters(parameters)
    
    modify = input("\nDo you want to modify any parameters? (Y/N): ")
    if modify.lower() == 'y':
        user_input = input("\nEnter parameters to modify (e.g., --seed 2 --nruns 3): ")
        parser = argparse.ArgumentParser()
        for key in parameters.keys():
            if key in ["popSizeVec", "inv_freq", "speciation", "demography", "invRange", "fixedSNPs", "tempRead", "nCarriers", "snpPositions"]:
                parser.add_argument(f"--{key}", nargs='+')  # Allow multiple values for these specific keys
            else:
                parser.add_argument(f"--{key}")  # Single value for other keys
          
        args = parser.parse_args(user_input.split())
        
        # Update the parameters dictionary with provided values
        for key, value in vars(args).items():
            if value is not None:
                if key in ["popSizeVec", "inv_freq", "speciation", "demography", "invRange", "fixedSNPs", "tempRead", "nCarriers", "snpPositions"]:
                    parameters[key] = " ".join(value)
                else:
                    parameters[key] = value

        print("\nUpdated Parameters:")
        print_parameters(parameters)
    else:
        print("\nNo changes made.")


# Running C++ program with parameters from Python
modify_params(parameters)
params=list(parameters.values())
run_cpp_program(params)

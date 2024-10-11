import subprocess

def run_cpp_program(params):
    # Call the compiled C++ program and pass the parameters
    command = ["./labp_v19"] + params
    result = subprocess.run( command, capture_output=True, text=True)

    # Print the output received from the C++ program
    print("C++ Program Output:\n", result.stdout)
    if result.stderr:
        print("Errors:\n", result.stderr)

## Default parameters
parameters = {
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
    modify = input("\nDo you want to modify any paramters? (Y/N): ")
    
    if modify == 'Y' or modify == 'y':
        while True:
            key = input("\nEnter the paramter name to modify (or type 'done' to finish): ")
            if key == 'done':
                break
            if key in parameters:
                new_value = input(f"Enter new value for {key} (current value: {parameters[key]}): ")
                parameters[key] = new_value
            else:
                print("Invalid parameter name. Please try again.")
        
        print("\nUpdated parameters:")
        print_parameters(parameters)
    else:
        print("\nNo changes made.")



# Running C++ program with parameters from Python
modify_params(parameters)
params=list(parameters.values())
run_cpp_program(params)

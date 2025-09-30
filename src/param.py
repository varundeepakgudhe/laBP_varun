import subprocess
import argparse

def run_cpp_program(params):
    # Call the compiled C++ program and pass the parameters
    command = ["./executables/labp_v31"] + params
    result = subprocess.run( command, capture_output=True, text=True)
    
    # Print the output received from the C++ program
    print("C++ Program Output:\n", result.stdout)
    if result.stderr:
        print(result.stderr)
        with open("Output_log.txt", "w") as log_file:
            log_file.write(result.stderr)

# Taking input in Python
# param1 = input("Enter first parameter: ")
# param2 = input("Enter second parameter: ")

# nruns = "1"
# kingman_coal = "1"
# drift_sim = "0"
# msOutput = "1"
# popSizeVec = "10000 10000"
# inv_freq = "0 0"
# speciation = "1 10000 0"
# demography = "0 0 0"
# inv_age = "0"
# migRate = "0"
# BasesPerMorgan = "1e8"
# randPhi = "0"
# phi = "0"
# invRange = "0 1e3"
# fixedSNPs = "0 1000 1e-6"
# n_SNPs = "1"
# snpPositions = "1000.01 1000.011"
# randomSample = "0"
# tempRead = "10 0"
# nCarriers = "10 0" #cross check these variable names

parameters = {
    "seed": "1",
    "nruns": "1",
    "kingman_coal": "1",
    "drift_sim": "0",
    "msOutput": "1",
    "popSizeVec": "10000 10000",
    "inv_freq": "0.2 0.4",
    "speciation": "1 10000 0.2", #the third index is ancestor freq.
    "demography": "1 10000 0",
    "inv_age": "0",
    "migRate": "0.02",
    "BasesPerMorgan": "1e8",
    "randPhi": "0",
    "phi": "0.2", ## should be greater than 0
    "invRange": "0 1e3",
    "fixedSNPs": "1 10 1e-6", ##here second index is nsNPS/nsites
    "randSNP": "0",
    "sites": "998 998.1 998.2 998.3 998.4 998.5 998.6 998.7 998.8 998.9",
    "randomSample": "0",
    "P1_samples": "10 0",
    "P2_samples": "10 0" #if you keep p1 as 5 5 and p2 as 5 5 and phi is 0, then the program never ends, coz they never migrate in population and they never stop the inverssion. its always going to be 5 5.
}

def check_sites(parameters):
    n_sites=int(parameters["fixedSNPs"].split(" ")[1])
    sites=len(parameters["sites"].split(" "))
    if parameters["randSNP"]=="0":
        if sites!=n_sites:
            raise ValueError(f"Expected {n_sites} site positions, but got {sites}.")
    elif parameters["randSNP"]=="1":
        if sites!=2:
            raise ValueError(f"Expected 2 site positions, but got {sites}.")
    else:
        raise ValueError(f"Expected value 0 or 1 for randSNP")



def print_parameters(parameters):
    for key, value in parameters.items():
        print(f"{key}: {value}")

# def modify_params(parameters):
#     print("Current Parameters:")
#     print_parameters(parameters)
#     modify = input("\nDo you want to modify any paramters? (Y/N): ")
    
#     if modify == 'Y' or modify == 'y':
#         while True:
#             key = input("\nEnter the paramter name to modify (or type 'done' to finish): ")
#             if key == 'done':
#                 break
#             if key in parameters:
#                 new_value = input(f"Enter new value for {key} (current value: {parameters[key]}): ")
#                 parameters[key] = new_value
#             else:
#                 print("Invalid parameter name. Please try again.")
        
#         print("\nUpdated parameters:")
#         print_parameters(parameters)
#     else:
#         print("\nNo changes made.")

def modify_params(parameters):
    print("Current Parameters:")
    print_parameters(parameters)
    
    modify = input("\nDo you want to modify any parameters? (Y/N): ")
    if modify.lower() == 'y':
        user_input = input("\nEnter parameters to modify (e.g., --seed 2 --nruns 3): ")
        parser = argparse.ArgumentParser()
        for key in parameters.keys():
            if key in ["popSizeVec", "inv_freq", "speciation", "demography", "invRange", "fixedSNPs", "P1_samples", "P2_samples", "sites"]:
                parser.add_argument(f"--{key}", nargs='+')  # Allow multiple values for these specific keys
            else:
                parser.add_argument(f"--{key}")  # Single value for other keys
          
        args = parser.parse_args(user_input.split())
        
        # Update the parameters dictionary with provided values
        for key, value in vars(args).items():
            if value is not None:
                if key in ["popSizeVec", "inv_freq", "speciation", "demography", "invRange", "fixedSNPs", "P1_samples", "P2_samples", "sites"]:
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
check_sites(parameters)
run_cpp_program(params)
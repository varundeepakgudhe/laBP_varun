import subprocess

def run_cpp_program(params):
    # Call the compiled C++ program and pass the parameters
    command = ["./labp_v19"] + params
    result = subprocess.run( command, capture_output=True, text=True)

    # Print the output received from the C++ program
    print("C++ Program Output:\n", result.stdout)
    if result.stderr:
        print("Errors:\n", result.stderr)


# Running C++ program with parameters from Python
params=["10000","1", "0", "1", "10000 10000", "0 0", "1 10000 0", "0 0 0", "0", "0", "1e8", "0", "0", "0 1e3", "0 1000 1e-6", "1", "1000.01 1000.011", "0", "10 0", "10 0"]
run_cpp_program(params)

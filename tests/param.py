import subprocess

def run_cpp_program(param1, param2):
    # Call the compiled C++ program and pass the parameters
    result = subprocess.run(["./labp_v18", param1, param2], capture_output=True, text=True) ## assuming running this from the folder tests 
    
    # Print the output received from the C++ program
    print("C++ Program Output:\n", result.stdout)
    if result.stderr:
        print("Errors:\n", result.stderr)


# Running C++ program with parameters from Python
run_cpp_program("1", "10000")

import subprocess

# Replace 'script.sh' with the path to your shell script
script_path = './test_solution.sh'

# Run the shell script
process = subprocess.Popen(['bash', script_path], stdout=subprocess.PIPE, text=True)

stdout_output, _ = process.communicate()

# Print or process the stdout output as needed
print(stdout_output.split("\n"))
stdout_output = stdout_output.split("\n")
print(f"SPEEDUP:  {int(stdout_output[-2])/int(stdout_output[-3])}")
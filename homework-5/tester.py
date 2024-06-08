import time
import subprocess


def run_command(command):
    try:
        # Use shell=True to run the command in the shell
        process = subprocess.Popen(
            command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, shell=True)

        # Wait for the process to complete and capture the output and error
        stdout, stderr = process.communicate()

        # Check if the process terminated with a non-zero exit code
        if process.returncode != 0:
            print(f"Command failed with return code {process.returncode}")
            return (stderr)
        else:
            # Print the captured output
            return (stdout)

    except Exception as e:
        print(f"An error occurred: {e}")


seed = 12345


seq_cmd = f"echo {seed} | ./sequential_implementation"
par_cmd = f"echo {seed} | ./student_submission"

start = time.time()
seq_res = run_command(seq_cmd)
mid = time.time()
par_res = run_command(par_cmd)
end = time.time()

if (seq_res != par_res):
    print(f"FAILED")
    print(seq_res)
    print(par_res)
else:
    print(f"PASSED")

print(start, mid, end)
print(f"SPEEDUP is (not reliable on my machine) {(mid-start)/(end-mid)}")

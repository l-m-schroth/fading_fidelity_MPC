"""
This file can be used for debugging purposes of the sweep
"""
#!/usr/bin/env python3
import os
import subprocess
from trunk.Euler_sweeps.trajectory_tracking_trunk_ee_closed_loop_euler import get_options_list

# Customize your results prefix and results folder here:
RESULTS_PREFIX = "results_16_link"
RESULTS_FOLDER = "results"

def main():
    # Call get_options_list() once to obtain the list of MPC options.
    options = get_options_list()
    num_options = len(options)
    print(f"Found {num_options} MPC options.")

    # Get the current working directory (assumed to be the parent folder of this script)
    cwd = os.getcwd()
    
    for i in range(num_options):
        # Build the command string to run the experiment locally.
        command = (
            f"python trajectory_tracking_trunk_ee_closed_loop_euler.py "
            f"--index {i} --results_prefix {RESULTS_PREFIX} --results_folder {RESULTS_FOLDER}"
        )
        print(f"Running experiment for MPC option {i} locally with command:\n{command}\n")
        # Run the command and wait for it to complete before moving on to the next.
        subprocess.run(command, shell=True)
        
if __name__ == "__main__":
    main()

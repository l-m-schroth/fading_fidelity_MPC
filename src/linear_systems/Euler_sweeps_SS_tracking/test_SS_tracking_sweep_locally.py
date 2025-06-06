"""
This file can be used to test the sweeping functionality locally. May require adaptation to your Euler setup.
"""

#!/usr/bin/env python3
import os
import subprocess
from linear_systems.Euler_sweeps_SS_tracking.SS_tracking_sweep_euler import get_options_list

RESULTS_FOLDER = "results"  # results subfolder inside current directory

def main():
    options = get_options_list()
    num_options = len(options)
    print(f"Found {num_options} experiment configurations.")

    # Get current working directory
    cwd = os.getcwd()

    # Ensure results folder exists
    os.makedirs(RESULTS_FOLDER, exist_ok=True)

    for i in range(num_options):
        # Define the command to execute
        command = [
            "python", f"{cwd}/SS_tracking_sweep_euler.py",
            "--index", str(i),
            "--results_folder", os.path.join(cwd, RESULTS_FOLDER)
        ]

        print(f"Running experiment {i} locally: {' '.join(command)}")
        
        # Run the experiment and capture output
        process = subprocess.run(command, capture_output=True, text=True)

        # Print output for debugging
        print(f"Experiment {i} completed with return code {process.returncode}")
        if process.stdout:
            print("STDOUT:", process.stdout)
        if process.stderr:
            print("STDERR:", process.stderr)

if __name__ == "__main__":
    main()

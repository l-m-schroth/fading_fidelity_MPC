"""
This script can be used to schedule the Euler sweep for constraint SS tracking on the Euler cluster, schedules one job per MPC.
May require some changes regarding your setup on Euler. Requires HPIPM to be installed.
"""

#!/usr/bin/env python3
import osR
import subprocess
from linear_systems.Euler_sweeps_SS_tracking.SS_tracking_sweep_euler import get_options_list
import time

RESULTS_FOLDER = "results"  # Adjust as needed

def main():
    options = get_options_list()
    num_options = len(options)
    print(f"Found {num_options} experiment configurations.")

    cwd = os.getcwd()  # path to Euler_sweeps_SS_tracking

    for i in range(num_options):
        job_name = f"ss_track_{i}"

        # Pause every 1000 iterations for euler
        if (i + 1) % 1000 == 0:
            time.sleep(15 * 60)  # Sleep for 15 minutes (300 seconds)

        # Build the SLURM script
        job_script = f"""#!/bin/bash
#SBATCH --mem-per-cpu=4G
#SBATCH --constraint=EPYC_7763
#SBATCH --job-name={job_name}
#SBATCH --tmp=4000
#SBATCH --output={job_name}.out
#SBATCH --error={job_name}.err
#SBATCH --time=1:00:00

# Source bashrc to load environment variables
source ~/.bashrc

# Source HPIPM environment variables
source ~/ASL_Masterthesis/Workspace/hpipm/examples/python/env.sh

# Activate virtual environment
source ~/ASL_Masterthesis/myenv/bin/activate

# Run the Python script
python {cwd}/SS_tracking_sweep_euler.py --index {i} --results_folder {cwd}/{RESULTS_FOLDER}
"""
        print(f"Submitting job for experiment option {i}...")
        subprocess.run(["sbatch"], input=job_script.encode("utf-8"))

if __name__ == "__main__":
    main()


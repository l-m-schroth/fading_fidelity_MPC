"""
This file can be used to initiate the Euler sweep for the computations of the unconstraint case
"""

#!/usr/bin/env python3
import os
import subprocess
from linear_systems.Euler_sweeps_LQR.Unconstrained_sweep_euler import get_options_list

RESULTS_FOLDER = "results"            # results subfolder inside Euler_sweeps

def main():
    options = get_options_list()
    num_options = len(options)
    print(f"Found {num_options} experiment configurations.")

    # Current working directory (Euler_sweeps)
    cwd = os.getcwd()

    for i in range(num_options):
        job_name = f"theory_sweep_{i}"

        # Build the SLURM script
        job_script = f"""#!/bin/bash
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name={job_name}
#SBATCH --tmp=4000
#SBATCH --output={job_name}.out
#SBATCH --error={job_name}.err
#SBATCH --time=5:00:00

source ~/ASL_Masterthesis/myenv/bin/activate

# Run the experiment
python {cwd}/Unconstrained_sweep_euler.py --index {i} --results_folder {cwd}/{RESULTS_FOLDER}

"""
        print(f"Submitting job for experiment option {i}...")
        subprocess.run(["sbatch"], input=job_script.encode("utf-8"))

if __name__ == "__main__":
    main()

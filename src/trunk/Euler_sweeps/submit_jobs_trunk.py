"""
This file can be used to initiate a sweep on the Euler cluster. We use the same CPU type and a single thread for a comparison as fair as possible
May need an adaption of the paths dependent on the repo setup on the cluster.
"""

#!/usr/bin/env python3
import os
import subprocess
from trunk.Euler_sweeps.trajectory_tracking_trunk_ee_closed_loop_euler import get_options_list

# results prefix and results folder:
RESULTS_PREFIX = "results_tracking_oval_ellipsoidal"
RESULTS_FOLDER = "results_19_03"

def main():
    # Call get_options_list() once to obtain the list of MPC options.
    options = get_options_list()
    num_options = len(options)
    print(f"Found {num_options} MPC options.")

    # Get the current working directory (assumed to be Euler_sweeps)
    cwd = os.getcwd()
    
    for i in range(num_options):
        # Use the job name as the working directory name.
        job_name = f"traj_track_mpc_{i}"
        work_dir = os.path.join(cwd, f"dir_{job_name}_{i}")
        # Build the SLURM job script as a multi-line string.
        job_script = f"""#!/bin/bash
#SBATCH --constraint=EPYC_7763
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name={job_name}
#SBATCH --tmp=4000

# Activate the virtual environment
source ~/ASL_Masterthesis/myenv/bin/activate

export ACADOS_SOURCE_DIR=~/ASL_Masterthesis/Workspace/acados
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/ASL_Masterthesis/Workspace/acados/lib

# Create and change to the unique working directory matching the job name
mkdir -p {work_dir}
cd {work_dir}
echo "Working directory: {work_dir}"

# Run the experiment for MPC option index {i} with custom prefix and results folder.
python {cwd}/trajectory_tracking_trunk_ee_closed_loop_euler.py --index {i} --results_prefix {RESULTS_PREFIX} --results_folder {RESULTS_FOLDER}

# Change back to the original directory and delete the working directory
cd {cwd}
rm -rf {work_dir}

# Deactivate the virtual environment
deactivate
"""
        print(f"Submitting job for MPC option {i}...")
        # Submit the job using subprocess; the job script is fed via stdin.
        subprocess.run(["sbatch"], input=job_script.encode('utf-8'))
        
if __name__ == "__main__":
    main()





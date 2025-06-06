#!/usr/bin/env python3
import argparse
import os
import pickle
import gc
from copy import deepcopy
import numpy as np
import casadi as ca
import mujoco

import os
import sys

# Add the parent directory to the Python path
parent_dir = os.path.dirname(os.path.dirname(__file__))
if parent_dir not in sys.path:
    sys.path.insert(0, parent_dir)

from trunk.trunk_mpc_2d_acados import TrunkMPC2DOptions, TrunkMPC2DEETracking
from trunk.Trajectory_finding_script import select_periodic_trajectory
from trunk.trunk_utils import get_ee_position, compute_q, compute_q_casadi
from trunk.plotting_utils_trajectory_tracking import plot_open_loop_plan_ee
from trunk.ODE_chains.ODE_utils import get_ode_params_dict
from trunk.trajectory_tracking_trunk_closed_loop import closed_loop_ee_tracking_mujoco, closed_loop_ee_tracking_acados
from trunk.trunk_utils import generate_exponential_step_sizes

def closed_loop_experiment_trunk_ee_euler(results_filename, ocp_options):
    """
    Run the closed-loop experiment for a given MPC option.
    The result is saved to results_filename.
    """
    # Experiment settings
    duration = 5.0        # seconds
    use_mujoco = False
    frequency = 1.0
    perturbation_variance = 0.0 # only works when mujoco is used

    # Load MuJoCo model
    parent_dir = os.path.dirname(os.path.dirname(__file__))
    xml_path = os.path.join(parent_dir, "models", "chain_models", "chain_16_links_expanded.xml")
    Mjmodel = mujoco.MjModel.from_xml_path(xml_path)
    Mjdata = mujoco.MjData(Mjmodel)

    # Generate trajectory using the first element of n from ocp_options
    trajectory_name = "trajectory_oval"
    n_val = ocp_options.n[0] 
    Phi_t_func, Phi_dot_t_func = select_periodic_trajectory(n=n_val, frequency=frequency, trajectory_name=trajectory_name, plot=False)
    T = 1 / frequency  
    time_steps = np.linspace(0, T, 100)
    trajectory = np.array([Phi_t_func(t_val) for t_val in time_steps])
    trajectory_x = trajectory[:, 0]
    trajectory_z = trajectory[:, 1]

    # Create the trunk MPC object for this experiment
    trunk_MPC = TrunkMPC2DEETracking(Phi_t_func, Phi_dot_t_func, ocp_options)

    try:
        if use_mujoco:
            _, ee_pos, costs, solve_times, SQP_iters, constraint_violations = closed_loop_ee_tracking_mujoco(
                n=ocp_options.n, n_high=ocp_options.n_high,
                Mjmodel=deepcopy(Mjmodel), Mjdata=deepcopy(Mjdata),
                duration=duration, framerate=None,
                create_video=False, plot_open_loop_plan_bool=False,
                trajectory_x=trajectory_x, trajectory_z=trajectory_z, trunkMPC=trunk_MPC, perturbation_variance=perturbation_variance
            )
        else:
            _, _, _, ee_pos, costs, solve_times, SQP_iters, constraint_violations = closed_loop_ee_tracking_acados(
                n=ocp_options.n, n_high=ocp_options.n_high,
                duration=duration, plot_open_loop_plan_bool=False,
                trajectory_x=trajectory_x, trajectory_z=trajectory_z, trunkMPC=trunk_MPC
            )
        results = {
            'options': ocp_options,
            'costs': costs,
            'solve_times': solve_times,
            'SQP_iters': SQP_iters,
            'constraint_violations': constraint_violations,
            'ee_pos': ee_pos
        }
    except Exception as e:
        print(f"Experiment failed for options {ocp_options}: {e}")
        results = None
    finally:
        del trunk_MPC
        gc.collect()

    if results is not None:
        with open(results_filename, 'wb') as f:
            pickle.dump(results, f)
        print(f"Results saved to {results_filename}")
    else:
        print("No results generated.")

def get_options_list():
    options_list = []
    dt_inital = 0.01
    exponential_rates = [1.0, 1.02, 1.04, 1.06, 1.08]
    QP_solvers = ['PARTIAL_CONDENSING_HPIPM', 'FULL_CONDENSING_HPIPM']
    n_values = [
        [16], [16, 8], [16, 4], [16, 8, 4], [16, 8, 'p'], [16, 8, 4, 'p'],    
    ]
    N_list_values = {
        "1": [[10], [15], [20], [25], [30], [35], [40], [45], [50]], # stages for standart ocp 
        "2": [[5, 5], [5, 10], [10, 5], [10, 10], [10, 20], [10, 30], [20, 20], [20, 30], [30, 20], [25, 25]], # stages for 2 phase problems
        "3": [[2, 5, 5], [5, 5, 5], [10, 5, 5], [5, 10, 10], [10, 10, 10], [10, 15, 15],
              [10, 15, 20], [20, 10 , 10]], # stages for 3-phase problems
        "4": [[5,5,5,5], [5, 5, 10, 10], [10, 10, 10, 10], [10, 10, 5, 5]] # stages for 4-phase problems
    }

    for nlp_solver_type in ["SQP"]:
        for n in n_values:
            for N_list in N_list_values[f"{len(n)}"]:
                for rate in exponential_rates:
                    for qp_solver in QP_solvers:
                        options_list.append(TrunkMPC2DOptions(
                            nlp_solver_type=nlp_solver_type,
                            n=n,
                            N_list=N_list,
                            n_high=16,
                            qp_solver = qp_solver,
                            dt=generate_exponential_step_sizes(initial_step=dt_inital, rate=rate, num_steps=sum(N_list), plot=False),
                            ub_x_ee=[0.2, 1e15]
                        ))

    return options_list

def run_experiment_for_index(index, results_folder, results_prefix):
    options = get_options_list()
    if index < 0 or index >= len(options):
        raise ValueError(f"Index {index} is out of range (0 to {len(options)-1}).")
    ocp_option = options[index]
    print(f"Running experiment for option index {index} with options: {ocp_option}")

    # Ensure the results folder exists.
    if not os.path.exists(results_folder):
        os.makedirs(results_folder)

    results_filename = os.path.join(results_folder, f"{results_prefix}_{index}.pkl")
    closed_loop_experiment_trunk_ee_euler(results_filename, ocp_option)

def main():
    parser = argparse.ArgumentParser(description="Run closed-loop experiment for the i-th MPC option.")
    parser.add_argument("--index", type=int, required=True, help="Index of the MPC option to run.")
    parser.add_argument("--results_prefix", type=str, default="results_tracking_oval_ellipses", help="File prefix for the saved results file.")
    parser.add_argument("--results_folder", type=str, default="results_19_03", help="Folder in which to save the results.")
    args = parser.parse_args()

    # Use the given results folder and prefix. If you want the folder to be inside Euler_sweeps,
    # you can convert it to an absolute path relative to this file.
    current_dir = os.path.dirname(__file__)
    results_folder = os.path.join(current_dir, args.results_folder)

    run_experiment_for_index(args.index, results_folder, args.results_prefix)

if __name__ == "__main__":
    main()


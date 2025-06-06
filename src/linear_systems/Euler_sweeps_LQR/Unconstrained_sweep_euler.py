#!/usr/bin/env python3
import os
import sys
import pickle
import argparse
import numpy as np
from scipy.linalg import solve_discrete_lyapunov

# If needed, insert the parent dir so we can import the "theory" folder or "utils" modules
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if parent_dir not in sys.path:
    sys.path.insert(0, parent_dir)

# Import everything you need from your project
from linear_systems.utils_linear_systems import get_LQR_info, compute_P_Lyapunov, get_number_optimization_vars
from linear_systems.utils_linear_systems import get_exponential_rates, get_model_schedules, divide_evenly
from linear_systems.utils_linear_systems import generate_augmented_linear_system, generate_exponential_step_sizes
import control as ct

# If you have your Q_tilde code in "theory" or "theory.utils_theory" – adjust as needed:
from linear_systems.utils_linear_systems import compute_Q_tilde_1_matrices_efficient, compute_Q_tilde

import pandas as pd

def get_options_list():
    """
    Create a base list of experiments (the same as before),
    but for each base experiment, generate multiple versions
    with different desired_lookahead values.
    """
    base_experiments = []

    # Example base experiment(s)
    # (Adapt to your liking, these are placeholders.)
    exp_base_1 = {
        "experiment_name": "n10_nu0",
        "random_seed": 42,
        "n": 10,
        "n_unstable": 0,
        "m": 4,
        "dt": 0.1,
        "model_order_min": 3,
        "N_steps_dt_increase_min": 4,
        "N_models_max": 3,
        "terminal_mode": "full",
    }
    base_experiments.append(exp_base_1)
    exp_base_2 = {
        "experiment_name": "n10_nu1",
        "random_seed": 42,
        "n": 10,
        "n_unstable": 1,
        "m": 4,
        "dt": 0.1,
        "model_order_min": 3,
        "N_steps_dt_increase_min": 3,
        "N_models_max": 3,
        "terminal_mode": "full",
    }
    base_experiments.append(exp_base_2)
    exp_base_3 = {
        "experiment_name": "n50_nu0",
        "random_seed": 42,
        "n": 50,
        "n_unstable": 0,
        "m": 10,
        "dt": 0.1,
        "model_order_min": 10,
        "N_steps_dt_increase_min": 4,
        "N_models_max": 3,
        "terminal_mode": "full",
    }
    base_experiments.append(exp_base_3)
    exp_base_4 = {
        "experiment_name": "n50_nu1",
        "random_seed": 42,
        "n": 50,
        "n_unstable": 1,
        "m": 10,
        "dt": 0.1,
        "model_order_min": 10,
        "N_steps_dt_increase_min": 4,
        "N_models_max": 3,
        "terminal_mode": "full",
    }
    base_experiments.append(exp_base_4)
    exp_base_5 = {
        "experiment_name": "n150_nu0",
        "random_seed": 42,
        "n": 150,
        "n_unstable": 0,
        "m": 15,
        "dt": 0.1,
        "model_order_min": 30,
        "N_steps_dt_increase_min": 4,
        "N_models_max": 3,
        "terminal_mode": "full",
    }
    base_experiments.append(exp_base_5)
    exp_base_6 = {
        "experiment_name": "n150_nu3",
        "random_seed": 42,
        "n": 150,
        "n_unstable": 3,
        "m": 15,
        "dt": 0.1,
        "model_order_min": 30,
        "N_steps_dt_increase_min": 4,
        "N_models_max": 3,
        "terminal_mode": "full",
    }
    base_experiments.append(exp_base_6)

    # The lookaheads we want to sweep over:
    lookaheads = np.linspace(0.5, 5.5, 20)

    all_options = []
    for base_exp in base_experiments:
        for lh in lookaheads:
            # Copy the base experiment dictionary
            exp = dict(base_exp)

            # Now set the lookahead
            exp["desired_lookahead"] = lh

            # Adjust experiment_name and/or prefix so it’s unique for each
            exp["experiment_name"] = f"{exp['experiment_name']}_lh{lh}"

            # Add to the list of full experiments
            all_options.append(exp)

    return all_options


def run_experiment(opts):
    """
    Perform the experiment logic for one set of options.
    This is where you replicate your snippet but use Q_tilde 
    for cost evaluation, etc.
    """

    # Unpack the dictionary
    random_seed = opts["random_seed"]
    n = opts["n"]
    n_unstable = opts["n_unstable"]
    m = opts["m"]
    desired_lookahead = opts["desired_lookahead"]
    dt = opts["dt"]
    model_order_min = opts["model_order_min"]
    N_models_max = opts["N_models_max"]
    terminal_mode = opts["terminal_mode"]
    A = opts["A"]
    B = opts["B"]

    # Set random seed
    np.random.seed(random_seed)

    # Generate the system
    # Note: l can be anything if you want to keep it as identity output
    # or pass in a param. We'll set l = n for a fully observed system.
    l = n
    # A, B, C, D = generate_augmented_linear_system(n, n_unstable, m, l,
    #                                              C_eye=True, discrete=False)
    # Ground truth discrete system for baseline LQR
    sys_gt_ct = ct.ss(A, B, np.eye(A.shape[0]), np.zeros((A.shape[0], B.shape[1])))

    # Q, R for LQR
    Q = np.eye(n)   # you might want to scale these
    R = 0.1*np.eye(m)

    # Get baseline "full-order" CARE solution
    P_full_CARE = ct.care(sys_gt_ct.A, sys_gt_ct.B, Q, R)[0]

    # Build the set of (model_orders, horizon_lengths, rates, etc.)
    N = int(desired_lookahead / dt)  # integer steps at dt
    N_steps_dt_increase_min = opts["N_steps_dt_increase_min"]
    exp_rates, exp_rate_steps = get_exponential_rates(N_steps_dt_increase_min, N, step=2)

    # Range of single-model orders
    model_orders = np.arange(model_order_min, n + 1)

    # "Mixed" model schedules
    schedules = get_model_schedules(n, list(model_orders),
                                    model_order_min, N_models_max, step=1)

    # For progress stats
    N_controllers_fixed_model = len(model_orders)*len(exp_rates)
    N_controllers_mixed = len(schedules) - len(model_orders)
    N_controllers_tot = N_controllers_fixed_model + N_controllers_mixed

    # We'll collect results here
    results_list = []

    progress_counter = 1
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # MAIN LOOP:
    #   1) For each schedule in "schedules"
    #      - If schedule has multiple chunks, we only do rate=1.0, 
    #        horizon_lengths = evenly divided
    #      - If single chunk, we do all exp_rates, horizon_lengths combos
    #   2) For each resulting setup, compute the LQR gains,
    #      then compute the "Q_tilde" cost for each x0
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Discretize the ground truth system 
    sys_gt_disc_step = ct.c2d(sys_gt_ct, dt, method='zoh')
    A_d = sys_gt_disc_step.A
    B_d = sys_gt_disc_step.B
    
    # Build Q_tilde      
    #Q1_1, Q1_2, Q1_3 = compute_Q_tilde_1_matrices_efficient(sys_gt_ct.A, sys_gt_ct.B, Q, dt)

    for schedule in schedules:
        num_chunks = len(schedule)
        # If more than one chunk -> "mixed model" => rate=1.0 only
        if num_chunks > 1:
            local_rates = [1.0]
            local_horizon_len_list = [divide_evenly(N, num_chunks)]
        else:
            # Single-model => test all exponential rates
            local_rates = exp_rates
            local_horizon_len_list = exp_rate_steps

        for rate, horizon_lens in zip(local_rates, local_horizon_len_list):
            # Build the step sizes
            step_sizes = generate_exponential_step_sizes(dt, rate, sum(horizon_lens))

            # 1) Build the LQR feedback from get_LQR_info
            K_list, _, _ = get_LQR_info(A, B, Q, R,
                                        schedule,
                                        horizon_lengths=horizon_lens,
                                        step_sizes=step_sizes,
                                        terminal_mode=terminal_mode)
            # For consistency, let's define K_applied as K_list[0].
            # In your code, you might re-check each step, but let's keep it simple.
            K_applied = K_list[0]
            
            # 2) compute Q_tilde for the chosen K (K_applied)
            #Q_tilde_mat = compute_Q_tilde(K_applied, Q1_1, Q1_2, Q1_3, R, dt)

            # Find corresponding P that satisfies DARE
            A_cl = A_d + B_d @ K_applied
            #P = solve_discrete_lyapunov(A_cl.T, Q_tilde_mat)
            P = compute_P_Lyapunov(A_d, B_d, K_applied, Q * dt, R * dt)

            # Number of variables (condensed / uncondensed)
            N_optim_vars_uncond = get_number_optimization_vars(schedule,
                                                               horizon_lens, m,
                                                               condesing=False)
            N_optim_vars_cond = get_number_optimization_vars(schedule,
                                                             horizon_lens, m,
                                                             condesing=True)

            results_list.append({
                "A":sys_gt_ct.A,
                "B":sys_gt_ct.B,
                "opts": opts,
                "schedule": schedule,
                "horizon_lengths": horizon_lens,
                "rate": rate,
                "N_optim_vars_uncondensed": N_optim_vars_uncond,
                "N_optim_vars_condensed": N_optim_vars_cond,
                "t_lookahead": desired_lookahead,
                "P_cont": P,
                "P_LQR": P_full_CARE,
                "is_stable": np.all(np.abs(np.linalg.eigvals(A_cl)) < 1)
            })

            print(f"[{opts['experiment_name']}] "
                  f"MPC {progress_counter}/{N_controllers_tot} done. "
                  f"(schedule={schedule}, rate={rate})")
            progress_counter += 1

    return results_list


def main():
    parser = argparse.ArgumentParser(
        description="Run a single multi-fidelity linear MPC experiment on Euler, identified by index."
    )
    parser.add_argument("--index", type=int, required=True,
                        help="Index of the experiment in the options list.")
    parser.add_argument("--results_folder", type=str, default="results",
                        help="Folder where results .pkl will be saved. Relative to current dir.")
    args = parser.parse_args()

    # Build the entire list of experiments
    all_opts = get_options_list()
    num_exps = len(all_opts)
    if args.index < 0 or args.index >= num_exps:
        raise ValueError(f"Invalid --index {args.index}. Only {num_exps} experiments available.")

    # Select that experiment
    exp_opts = all_opts[args.index]

    # Run it
    results_list = run_experiment(exp_opts)

    # Save results
    filename = f"{exp_opts['experiment_name']}_results_{args.index}.pkl"
    results_path = os.path.join(args.results_folder, filename)

    os.makedirs(args.results_folder, exist_ok=True)

    with open(results_path, "wb") as f:
        pickle.dump(results_list, f)

    print(f"Saved results to {results_path}")


if __name__ == "__main__":
    main()

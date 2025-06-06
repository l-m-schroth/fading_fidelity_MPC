#!/usr/bin/env python3
import os
import sys
import pickle
import argparse
import numpy as np
import control as ctrl

# Insert parent dir for your local imports
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if parent_dir not in sys.path:
    sys.path.insert(0, parent_dir)

# Local/project imports
from linear_systems.utils_linear_systems import generate_augmented_linear_system
from linear_systems.utils_linear_systems import get_exponential_rates, get_model_schedules, divide_evenly
from linear_systems.utils_linear_systems import compute_P_Lyapunov
from linear_systems.LinearSSTrackingMPC_HPIPM import LinearSSTrackingMPChpipm, LinearSSTrackingMPCOptions
from linear_systems.LinearSSTrackingMPC_HPIPM_condensed import LinearSSTrackingMPChpipmCondensing
from linear_systems.utils_linear_systems import closed_loop_simulation_2  
from scipy.linalg import solve_discrete_are

def get_options_list():
    """
    Builds and returns a list of *all* experiment configurations. 
    Each entry in the returned list is one job -> one unique MPC run.

    """
                                      
    base_experiments = []

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Example: define your base experiments here
    # (Adjust or add more as you wish)
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # exp_base_1 = {
    #     "experiment_name": "n50_nu0_nonoise",
    #     "random_seed": 47,
    #     "n": 50,
    #     "n_unstable": 0,
    #     "m": 5,
    #     "dt": 0.05,
    #     "model_order_min": 10,
    #     "N_steps_dt_increase_min": 4,
    #     "N_models_max": 3,
    #     # Per-experiment constraint:
    #     "A_constr": np.hstack([np.eye(3), np.zeros((3,47))]),
    #     "b_constr": np.array([2.2]*3, dtype=float),
    #     # The initial state
    #     "x0": np.ones(50),
    #     # The total simulation duration in seconds
    #     "duration": 15.0,  
    #     "error_scale": 0.0,
    # }
    # base_experiments.append(exp_base_1)
    exp_base_2 = {
        "experiment_name": "n50_nu1_nonoise",
        "random_seed": 41,
        "n": 50,
        "n_unstable": 1,
        "m": 5,
        "dt": 0.05,
        "model_order_min": 10,
        "N_steps_dt_increase_min": 4,
        "N_models_max": 3,
        # Per-experiment constraint:
        "A_constr": np.hstack([np.eye(3), np.zeros((3,47))]),
        "b_constr": np.array([7]*3, dtype=float),
        # The initial state
        "x0": np.ones(50),
        # The total simulation duration in seconds
        "duration": 20.0,
        "error_scale": 0.0,
    }
    base_experiments.append(exp_base_2)
    # dt_next = 0.05
    # exp_base_3 = {
    #     "experiment_name": "n50_nu0_noise_dt0.05",
    #     "random_seed": 47,
    #     "n": 50,
    #     "n_unstable": 0,
    #     "m": 5,
    #     "dt": dt_next,
    #     "model_order_min": 10,
    #     "N_steps_dt_increase_min": 4,
    #     "N_models_max": 3,
    #     # Per-experiment constraint:
    #     "A_constr": np.hstack([np.eye(3), np.zeros((3,47))]),
    #     "b_constr": np.array([2.2]*3, dtype=float),
    #     # The initial state
    #     "x0": np.ones(50),
    #     # The total simulation duration in seconds
    #     "duration": 15.0,
    #     "error_scale": dt_next*2 
    # }
    # base_experiments.append(exp_base_3)
    # exp_base_4 = {
    #     "experiment_name": "n50_nu1_noise_dt0.05",
    #     "random_seed": 41,
    #     "n": 50,
    #     "n_unstable": 1,
    #     "m": 5,
    #     "dt": dt_next,
    #     "model_order_min": 10,
    #     "N_steps_dt_increase_min": 4,
    #     "N_models_max": 3,
    #     # Per-experiment constraint:
    #     "A_constr": np.hstack([np.eye(3), np.zeros((3,47))]),
    #     "b_constr": np.array([10]*3, dtype=float),
    #     # The initial state
    #     "x0": np.ones(50),
    #     # The total simulation duration in seconds
    #     "duration": 15.0,  
    #     "error_scale": dt_next*2 
    # }
    # base_experiments.append(exp_base_4)
    # dt_next = 0.1
    # exp_base_5 = {
    #     "experiment_name": "n50_nu0_noise_dt0.1",
    #     "random_seed": 47,
    #     "n": 50,
    #     "n_unstable": 0,
    #     "m": 5,
    #     "dt": dt_next,
    #     "model_order_min": 10,
    #     "N_steps_dt_increase_min": 4,
    #     "N_models_max": 3,
    #     # Per-experiment constraint:
    #     "A_constr": np.hstack([np.eye(3), np.zeros((3,47))]),
    #     "b_constr": np.array([2.2]*3, dtype=float),
    #     # The initial state
    #     "x0": np.ones(50),
    #     # The total simulation duration in seconds
    #     "duration": 15.0,
    #     "error_scale": dt_next*2 
    # }
    # base_experiments.append(exp_base_5)
    # exp_base_6 = {
    #     "experiment_name": "n50_nu1_noise_dt0.1",
    #     "random_seed": 41,
    #     "n": 50,
    #     "n_unstable": 1,
    #     "m": 5,
    #     "dt": dt_next,
    #     "model_order_min": 10,
    #     "N_steps_dt_increase_min": 4,
    #     "N_models_max": 3,
    #     # Per-experiment constraint:
    #     "A_constr": np.hstack([np.eye(3), np.zeros((3,47))]),
    #     "b_constr": np.array([10]*3, dtype=float),
    #     # The initial state
    #     "x0": np.ones(50),
    #     # The total simulation duration in seconds
    #     "duration": 15.0,  
    #     "error_scale": dt_next*2 
    # }
    # base_experiments.append(exp_base_6)
    # dt_next = 0.01
    # exp_base_7 = {
    #     "experiment_name": "n50_nu0_noise_dt0.01",
    #     "random_seed": 47,
    #     "n": 50,
    #     "n_unstable": 0,
    #     "m": 5,
    #     "dt": dt_next,
    #     "model_order_min": 10,
    #     "N_steps_dt_increase_min": 4,
    #     "N_models_max": 3,
    #     # Per-experiment constraint:
    #     "A_constr": np.hstack([np.eye(3), np.zeros((3,47))]),
    #     "b_constr": np.array([2.2]*3, dtype=float),
    #     # The initial state
    #     "x0": np.ones(50),
    #     # The total simulation duration in seconds
    #     "duration": 15.0,
    #     "error_scale": dt_next*2 
    # }
    # base_experiments.append(exp_base_7)
    # exp_base_8 = {
    #     "experiment_name": "n50_nu1_noise_dt0.01",
    #     "random_seed": 41,
    #     "n": 50,
    #     "n_unstable": 1,
    #     "m": 5,
    #     "dt": dt_next,
    #     "model_order_min": 10,
    #     "N_steps_dt_increase_min": 4,
    #     "N_models_max": 3,
    #     # Per-experiment constraint:
    #     "A_constr": np.hstack([np.eye(3), np.zeros((3,47))]),
    #     "b_constr": np.array([10]*3, dtype=float),
    #     # The initial state
    #     "x0": np.ones(50),
    #     # The total simulation duration in seconds
    #     "duration": 15.0,  
    #     "error_scale": dt_next*2 
    # }
    # base_experiments.append(exp_base_8)
    # dt_next = 0.025
    # exp_base_9 = {
    #     "experiment_name": "n50_nu0_noise_dt0.025",
    #     "random_seed": 47,
    #     "n": 50,
    #     "n_unstable": 0,
    #     "m": 5,
    #     "dt": dt_next,
    #     "model_order_min": 10,
    #     "N_steps_dt_increase_min": 4,
    #     "N_models_max": 3,
    #     # Per-experiment constraint:
    #     "A_constr": np.hstack([np.eye(3), np.zeros((3,47))]),
    #     "b_constr": np.array([2.2]*3, dtype=float),
    #     # The initial state
    #     "x0": np.ones(50),
    #     # The total simulation duration in seconds
    #     "duration": 15.0,
    #     "error_scale": dt_next*2 
    # }
    # base_experiments.append(exp_base_9)
    # exp_base_10 = {
    #     "experiment_name": "n50_nu1_noise_dt0.025",
    #     "random_seed": 41,
    #     "n": 50,
    #     "n_unstable": 1,
    #     "m": 5,
    #     "dt": dt_next,
    #     "model_order_min": 10,
    #     "N_steps_dt_increase_min": 4,
    #     "N_models_max": 3,
    #     # Per-experiment constraint:
    #     "A_constr": np.hstack([np.eye(3), np.zeros((3,47))]),
    #     "b_constr": np.array([10]*3, dtype=float),
    #     # The initial state
    #     "x0": np.ones(50),
    #     # The total simulation duration in seconds
    #     "duration": 15.0,  
    #     "error_scale": dt_next*2 
    # }
    # base_experiments.append(exp_base_10)
    # If you have more base experiments, add them similarly:
    # exp_base_2 = { ... }
    # base_experiments.append(exp_base_2)
    # etc.

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # For each base experiment, we sweep over a set of lookaheads
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    lookaheads = [1.0, 1.5, 2.0, 2.5]

    all_jobs = []
    for base_exp in base_experiments:
        dt = base_exp["dt"]

        for lh in lookaheads:
            # 1) Copy the base experiment
            exp_cfg = dict(base_exp)
            exp_cfg["lookahead_distance"] = lh

            # 2) Compute integer horizon N from dt and lookahead
            N = int(np.ceil(lh / dt))

            # 3) Generate exponential rates + steps
            exp_rates, exp_rate_steps = get_exponential_rates(
                base_exp["N_steps_dt_increase_min"], 
                N, 
                step=2
            )

            # 4) Generate model schedules
            model_orders_fixed = np.arange(base_exp["model_order_min"], base_exp["n"] + 1)
            schedules = get_model_schedules(
                max_order_init=base_exp["n"],
                model_orders_fixed_model=model_orders_fixed.tolist(),
                absolute_min_model_order=base_exp["model_order_min"],
                N_models_max=base_exp["N_models_max"],
                step=2
            )

            # 5) For each schedule, figure out step sizes from the rates
            for schedule in schedules:
                num_chunks = len(schedule)
                if num_chunks > 1:
                    # "mixed" => single rate=1.0
                    local_rates = [1.0]
                    local_horizon_len_list = [divide_evenly(N, num_chunks)]
                else:
                    # single-model => test all exponential rates
                    local_rates = exp_rates
                    local_horizon_len_list = exp_rate_steps

                for rate, horizon_lens in zip(local_rates, local_horizon_len_list):
                    # Build an HPC job config with uncondensed HPIPM
                    from linear_systems.utils_linear_systems import generate_exponential_step_sizes
                    step_sizes = generate_exponential_step_sizes(dt, rate, sum(horizon_lens))

                    job_cfg = dict(exp_cfg)
                    job_cfg["schedule"] = schedule
                    job_cfg["horizon_lengths"] = horizon_lens
                    job_cfg["step_sizes"] = step_sizes
                    job_cfg["rate"] = rate
                    job_cfg["is_condensed"] = False  # default

                    # Unique job name
                    job_cfg["experiment_name"] = (
                        f"{base_exp['experiment_name']}_lh{lh}"
                        f"_sched{schedule}_rate{rate}_condFalse"
                    )
                    all_jobs.append(job_cfg)

                    # 6) If this schedule is the "exact full-order model" schedule, 
                    #    also create a second job for the condensed solver
                    #    i.e. single chunk of size n == base_exp["n"]
                    if (num_chunks == 1) and (schedule[0] == base_exp["n"]):
                        # This is the full-order schedule => add condensed version
                        cond_job_cfg = dict(job_cfg)
                        cond_job_cfg["is_condensed"] = True
                        cond_job_cfg["experiment_name"] = (
                            f"{base_exp['experiment_name']}_lh{lh}"
                            f"_sched{schedule}_rate{rate}_condTrue"
                        )
                        all_jobs.append(cond_job_cfg)

    return all_jobs


def run_experiment(opts):
    """
    Runs one closed-loop simulation of HPIPM-based SS Tracking MPC 
    (condensed or uncondensed) using the specified schedule, horizon, etc.
    """
    # Unpack
    random_seed = opts["random_seed"]
    n = opts["n"]
    n_unstable = opts["n_unstable"]
    m = opts["m"]
    dt = opts["dt"]
    x0 = opts["x0"]
    duration = opts["duration"]
    schedule = opts["schedule"]
    horizon_lengths = opts["horizon_lengths"]
    step_sizes = opts["step_sizes"]

    np.random.seed(random_seed)

    # Ground truth system
    A_ct, B_ct, C_ct, D_ct = generate_augmented_linear_system(
        n, n_unstable, m, l=n, C_eye=True, discrete=False
    )
    sys_ct = ctrl.ss(A_ct, B_ct, C_ct, D_ct)

    # Q, R
    Q = 10*np.eye(n)
    R = 0.01*np.eye(m)

    # Terminal cost from DARE
    sys_d = ctrl.c2d(sys_ct, dt, method='zoh')
    A_d = sys_d.A
    B_d = sys_d.B
    P_terminal = solve_discrete_are(A_d, B_d, Q * dt, R * dt)

    # Build MPC options
    from linear_systems.LinearSSTrackingMPC import LinearSSTrackingMPCOptions
    mpc_opts = LinearSSTrackingMPCOptions(
        model_orders=schedule,
        horizon_lengths=horizon_lengths,
        step_sizes=step_sizes,
        Q=Q,
        R=R,
        K=P_terminal,
        ground_truth_ct_system=sys_ct,
        A_constr=opts["A_constr"],
        b_constr=opts["b_constr"],
        soft_constraints=False,
        lam=0
    )

    # Create solver (condensed or not)
    if opts["is_condensed"]:
        from linear_systems.LinearSSTrackingMPC_HPIPM_condensed import LinearSSTrackingMPChpipmCondensing
        mpc = LinearSSTrackingMPChpipmCondensing(mpc_opts)
    else:
        mpc = LinearSSTrackingMPChpipm(mpc_opts)

    # Run closed-loop simulation
    simulation_freq = 1.0 / dt
    control_freq = 1.0 / dt
    closed_loop_data = closed_loop_simulation_2(
        simulation_freq,
        control_freq,
        mpc,
        x0,
        duration,
        error_scale=opts["error_scale"]
    )

    # Solve times
    comp_times = closed_loop_data["computation_times"]
    mean_solve_time = float(np.mean(comp_times)) if comp_times else 0.0
    std_solve_time = float(np.std(comp_times)) if comp_times else 0.0

    # Summation of cost
    total_cost = float(np.sum(closed_loop_data["costs_total"]))

    results = {
        "experiment_opts": opts,         # The entire experiment config
        "A_ground_truth": A_ct,
        "B_ground_truth": B_ct,
        "closed_loop_trajectory": closed_loop_data,
        "total_closed_loop_cost": total_cost,
        "mean_solve_time": mean_solve_time,
        "std_solve_time": std_solve_time,
        "is_condensed": opts["is_condensed"]
    }

    return results


def main():
    parser = argparse.ArgumentParser(
        description="Run a single SS tracking experiment (HPIPM-based) identified by index."
    )
    parser.add_argument("--index", type=int, required=True,
                        help="Index of the experiment in the options list.")
    parser.add_argument("--results_folder", type=str, default="results",
                        help="Folder where results .pkl will be saved.")
    args = parser.parse_args()

    # Build all possible experiment configurations
    all_configs = get_options_list()
    num_exps = len(all_configs)

    if args.index < 0 or args.index >= num_exps:
        raise ValueError(f"Invalid --index {args.index}. Only {num_exps} experiments available.")

    exp_opts = all_configs[args.index]

    # Run it
    result_data = run_experiment(exp_opts)

    # Save results
    # The user wants a filename that includes the lookahead distance, schedule, cond/no-cond, etc.
    # We'll keep it simple: we already have "experiment_name" in the config.
    job_name = exp_opts["experiment_name"]
    filename = f"{job_name}_results_{args.index}.pkl"

    os.makedirs(args.results_folder, exist_ok=True)
    save_path = os.path.join(args.results_folder, filename)

    with open(save_path, "wb") as f:
        pickle.dump(result_data, f)

    print(f"[INDEX {args.index}] Saved results to: {save_path}")


if __name__ == "__main__":
    main()

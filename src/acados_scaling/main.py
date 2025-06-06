"""
This file can be used to initiate a sweep over acados solver runs with over different parameter values and randomly generated linear systems.
We use it the plot the solver timing with the number of states n_x and the horizon length N (Figure 2.1).
"""

from acados_template import AcadosOcp, AcadosOcpSolver, AcadosModel
from acados_scaling.model import export_model
import numpy as np
import casadi as ca
import argparse
import yaml
from itertools import product
import pandas as pd
import os
import gc
import shutil
from datetime import datetime
from utils_shared import get_dir

STATS = ['time_tot',
        'time_lin',
        'time_sim',
        'time_sim_ad',
        'time_sim_la',
        'time_qp',
        'time_qp_solver_call',
        'time_qp_xcond',
        'time_glob',
        'time_solution_sensitivities',
        'time_reg',
        'time_preparation',
        'time_feedbacks'
] + [
        'sqp_iter',
        'ddp_iter',
        'nlp_iter',
        'qp_stat',
        'qp_iter',
        'statistics',
        'stat_m',
        'stat_n',
        'residuals',
        'alpha',
        'res_eq_all',
        'res_stat_all',
    ]


def load_sweep_config(file_path):
    with open(file_path, "r") as file:
        return yaml.safe_load(file)
    

def main(args):
    ocp = AcadosOcp()

    model = export_model(args.nx, args.nu)
    ocp.model = model

    nx = model.x.rows()
    nu = model.u.rows()
    N = args.N

    ocp.solver_options.N_horizon = N
    ocp.solver_options.tf = 1 # irrelevant due to discrete dynamics expression

    Q_mat = np.diag([1e1]*nx)
    R_mat = np.diag([1e-4]*nu)

    ocp.cost.cost_type = 'NONLINEAR_LS'
    ocp.model.cost_y_expr = ca.vertcat(model.x, model.u)
    ocp.cost.yref = np.zeros((nx+nu,))
    ocp.cost.W = ca.diagcat(Q_mat, R_mat).full()

    ocp.cost.cost_type_e = 'NONLINEAR_LS'
    ocp.cost.yref_e = np.zeros((nx,))
    ocp.model.cost_y_expr_e = model.x
    ocp.cost.W_e = Q_mat

    # Fmax = 80
    # ocp.constraints.lbu = np.array([-Fmax]*nu)
    # ocp.constraints.ubu = np.array([+Fmax]*nu)
    # ocp.constraints.idxbu = np.arange(nu)
    ocp.constraints.x0 = np.array([1.0]*nx)

    ocp.solver_options.qp_solver = args.qp_solver
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'
    ocp.solver_options.integrator_type = 'DISCRETE'
    ocp.solver_options.nlp_solver_type = 'SQP'

    ocp_solver = AcadosOcpSolver(ocp, generate=True, build=True)

    stats_list = []
    for _ in range(50): 
        stats = solve(ocp_solver, nx, args)
        stats_list.append(stats)

    del ocp_solver
    gc.collect()
    shutil.rmtree("c_generated_code", ignore_errors=True)

    return stats_list


def solve(ocp_solver, nx, args=None):
    x0 = np.random.rand(nx)
    ocp_solver.set(0, "lbx", x0)
    ocp_solver.set(0, "ubx", x0)
    status = ocp_solver.solve()
    ocp_solver.print_statistics()

    stats_dict = {}
    for stat in STATS:
        try:
            stats_dict[stat] = ocp_solver.get_stats(stat)
        except Exception:
            pass

    if status != 0:
        raise Exception(f"acados returned status {status}.")

    # Attach parameters used for this run
    if args:
        stats_dict.update({
            "N": args.N,
            "nx": args.nx,
            "nu": args.nu,
            "qp_solver": args.qp_solver,
        })

    return stats_dict


def sweep(args, sweep_config):
    results = []
    param_grid = sweep_config["parameters"]
    N_vals = param_grid["N"]["values"]
    nx_vals = param_grid["nx"]["values"]
    nu_vals = param_grid["nu"]["values"]
    qp_solver_vals = param_grid["qp_solver"]["values"]

    for N, nx, nu, qp_solver in product(N_vals, nx_vals, nu_vals, qp_solver_vals):
        args.N = N
        args.nx = nx
        args.nu = nu
        args.qp_solver = qp_solver

        print(f"Running sweep with N={N}, nx={nx}, nu={nu}, qp_solver={qp_solver}")
        run_results = main(args)
        time_tot_array = [r["time_tot"] for r in run_results] 
        results.append({
            "N": N,
            "nx": nx,
            "nu": nu,
            "qp_solver": qp_solver,
            "time_tot": time_tot_array
        })

    df = pd.DataFrame(results)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    data_dir = get_dir("data/acados_scaling")
    output_file = data_dir / f"sweep_results_{timestamp}.csv"
    df.to_csv(output_file, index=False)
    print(f"Sweep results saved to {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--N", type=int, default=20, help="Prediction horizon")
    parser.add_argument("--nx", type=int, default=10, help="Number of states")
    parser.add_argument("--nu", type=int, default=10, help="Number of controls")
    parser.add_argument("--qp-solver", type=str, default="FULL_CONDENSING_HPIPM", help="QP solver")
    parser.add_argument("--sweep", action="store_true", help="Run sweep")

    args = parser.parse_args()

    if args.sweep:
        sweep_config = load_sweep_config("sweep_config.yaml")
        sweep(args, sweep_config)
    else:
        main(args)


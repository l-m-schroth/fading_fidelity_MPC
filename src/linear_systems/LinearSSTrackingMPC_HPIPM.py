import os
import numpy as np
import control as ct
from hpipm_python import *
from hpipm_python.common import *
from dataclasses import dataclass
import time

# Import the base class
from linear_systems.mor.balancedTruncationUnstable import balanced_truncation_unstable
from linear_systems.LinearSSTrackingMPC import LinearSSTrackingMPC, LinearSSTrackingMPCOptions

class LinearSSTrackingMPChpipm(LinearSSTrackingMPC):
    """
    Linear State-Space Tracking MPC using HPIPM as the solver backend.
    Inherits from the CVXPY-based implementation to reuse model reduction functionality,
    but overrides the solve() method to use HPIPM.
    Take in LinearSSTrackingMPCOptions, where the model and step size schedule can be specified.
    """
    def __init__(self, opts: LinearSSTrackingMPCOptions):
        super().__init__(opts)

    def build_MPC_problem(self):
        N_tot = sum(self.opts.horizon_lengths)
        
        # 1) Create OCP QP dimension object for an N_tot-stage problem
        self.dim = hpipm_ocp_qp_dim(N_tot)

        # 2) Set dimensions stage by stage
        for k in range(N_tot + 1):
            if k < N_tot:
                i = self.chunk_index(k)
                n_red_i = self.opts.model_orders[i]
                self.dim.set('nx', n_red_i, k)
            else:
                # Terminal stage
                self.dim.set('nx', self.opts.model_orders[-1], N_tot)

            if k < N_tot:
                self.dim.set('nu', self.m, k)
            
            # Handle constraints
            if self.opts.A_constr is not None and self.opts.b_constr is not None:
                n_constr = self.opts.b_constr.shape[0]
                # If you want them as general constraints:
                if self.opts.soft_constraints:
                    self.dim.set('nsg', n_constr, k)
                    self.dim.set('ng', n_constr, k)
                    self.dim.set('ns', n_constr, k)
                else:
                    self.dim.set('ng', n_constr, k)

        # -----------------------------------------------------
        # >>> Enforce initial state at stage 0 as box constraints <<<
        # Make sure all states at stage 0 are fixed to x0_red
        self.dim.set('nbx', self.opts.model_orders[0], 0)
        # -----------------------------------------------------

        # print to shell
        self.dim.print_C_struct()

        # 3) Build QP
        self.qp = hpipm_ocp_qp(self.dim)

        # >>> Instead of self.qp.set('Jx', ...), do:
        # Identity for the entire reduced state
        Jx0 = np.eye(self.opts.model_orders[0])
        self.qp.set('Jbx', Jx0, 0)

        # The actual numeric bounds for x0 will be set in solve(), once x0 is known
        # Initialize them to zero for now
        self.qp.set('lbx', np.zeros((self.opts.model_orders[0],)), 0)
        self.qp.set('ubx', np.zeros((self.opts.model_orders[0],)), 0)
        # -----------------------------------------------------

        # Fill in the rest (dynamics, costs, constraints) for k in 0..N_tot-1
        for k in range(N_tot):
            i = self.chunk_index(k)
            n_red_i = self.opts.model_orders[i]
            dt_k = self.opts.step_sizes[k]

            # Retrieve discrete-time reduced model
            A_tilde, B_tilde, _, W, V, _ = self.reduced_models[(n_red_i, dt_k)]

            # Bridging for next step
            if k < N_tot - 1:
                j = self.chunk_index(k+1)
                n_red_j = self.opts.model_orders[j]
                dt_kplus = self.opts.step_sizes[k+1]
                _, _, _, W_next, V_next, _ = self.reduced_models[(n_red_j, dt_kplus)]
                Psi_ij = W_next.T @ V   # bridging
            else:
                Psi_ij = np.eye(self.opts.model_orders[-1])

            # Set state-update matrices
            self.qp.set('A', Psi_ij @ A_tilde, k)
            self.qp.set('B', Psi_ij @ B_tilde, k)
            self.qp.set('b', np.zeros((Psi_ij.shape[0],)), k)

            # Stage cost
            Q_scaled = dt_k * (V.T @ self.opts.Q @ V)
            R_scaled = dt_k * self.opts.R
            self.qp.set('Q', Q_scaled, k)
            self.qp.set('R', R_scaled, k)
            self.qp.set('q', np.zeros(n_red_i), k)
            self.qp.set('r', np.zeros(self.m), k)

            # Constraints
            if self.opts.A_constr is not None and self.opts.b_constr is not None:
                A_constr_red = self.opts.A_constr @ V
                n_constr = self.opts.b_constr.shape[0]
                if self.opts.soft_constraints:
                    #idxs = np.arange(n_constr)  # Assuming all constraints are softened
                    #self.qp.set('idxs', idxs, k)
                    # Soft constraints => 'C', 'lg', 'ug' + cost on slack
                    self.qp.set('C', A_constr_red, k)
                    self.qp.set('D', np.zeros((n_constr, self.m)), k)
                    self.qp.set('lg', -self.opts.b_constr.reshape(n_constr,1), k)
                    self.qp.set('ug', self.opts.b_constr.reshape(n_constr,1), k)

                    # Slack cost
                    Zl = np.diag(np.zeros(n_constr).reshape(n_constr,1))#np.diag((dt_k * self.opts.lam * np.ones(n_constr)).reshape(n_constr,1))
                    Zu = np.diag(np.zeros(n_constr).reshape(n_constr,1))#np.diag((dt_k * self.opts.lam * np.ones(n_constr)).reshape(n_constr,1))
                    zl = np.zeros((n_constr, 1)) #(dt_k * self.opts.lam * np.ones(n_constr))
                    zu = (dt_k * self.opts.lam * np.ones(n_constr))#.reshape(n_constr,1) # dt_k * self.opts.lam * np.ones(n_constr) np.zeros(n_constr).reshape(n_constr,1)#(dt_k * self.opts.lam * np.ones(n_constr)).reshape(n_constr,1) # dt_k * self.opts.lam * np.ones(n_constr)
                    Jsg = np.eye(n_constr)

                    self.qp.set('Zl', Zl, k)
                    self.qp.set('Zu', Zu, k)
                    self.qp.set('zl', zl, k)
                    self.qp.set('zu', zu, k)
                    self.qp.set('Jsg', Jsg, k)

                    slack_bounds = np.zeros((n_constr,1))
                    self.qp.set('lus', slack_bounds, k)
                    self.qp.set('lls', slack_bounds, k)
                else:
                    # Hard constraints => 'C', 'lg', 'ug' 
                    self.qp.set('C', A_constr_red, k)
                    self.qp.set('D', np.zeros((n_constr, self.m)), k)
                    self.qp.set('lg', -1000*self.opts.b_constr.reshape(n_constr,1), k)
                    self.qp.set('ug', self.opts.b_constr.reshape(n_constr,1), k)

                    Jsg = np.eye(n_constr)
                    self.qp.set('Jsg', Jsg, k)

        # Terminal cost
        dt_end = self.opts.step_sizes[-1]
        n_red_end = self.opts.model_orders[-1]
        _, _, _, _, V_end, _ = self.reduced_models[(n_red_end, dt_end)]
        Q_terminal = V_end.T @ self.opts.K @ V_end
        self.qp.set('Q', Q_terminal, N_tot)
        self.qp.set('q', np.zeros(n_red_end), N_tot)

        # print to shell
        self.qp.print_C_struct()

        # Create QP solution object
        self.qp_sol = hpipm_ocp_qp_sol(self.dim)

        # Solver arguments
        self.arg = hpipm_ocp_qp_solver_arg(self.dim, 'robust')
        self.arg.set('mu0', 1e4)
        self.arg.set('iter_max', 200)
        self.arg.set('tol_stat', 1e-6)
        self.arg.set('tol_eq', 1e-6)
        self.arg.set('tol_ineq', 1e-6)
        self.arg.set('tol_comp', 1e-6)
        self.arg.set('reg_prim', 1e-12)

        # Initialize solver
        self.solver = hpipm_ocp_qp_solver(self.dim, self.arg)

    def solve(self, x0: np.ndarray):
        """
        Solve the MPC problem with the given initial state using HPIPM.
        """
        # Store the initial full state
        self.x0 = x0.copy()

        # Compute reduced initial state
        dt_0 = self.opts.step_sizes[0]
        n_red_0 = self.opts.model_orders[0]
        _, _, _, W_0, V_0, _ = self.reduced_models[(n_red_0, dt_0)]
        x0_reduced = W_0.T @ x0

        # Update the initial-state box constraints for stage 0
        self.qp.set('lbx', x0_reduced, 0)
        self.qp.set('ubx', x0_reduced, 0)

        # Solve
        start_time = time.time()
        self.solver.solve(self.qp, self.qp_sol)
        solve_time = time.time() - start_time

        # Retrieve and print residuals
        res_stat = self.solver.get('max_res_stat')
        res_eq = self.solver.get('max_res_eq')
        res_ineq = self.solver.get('max_res_ineq')
        res_comp = self.solver.get('max_res_comp')
        print(f"Residuals - res_stat: {res_stat}, res_eq: {res_eq}, res_ineq: {res_ineq}, res_comp: {res_comp}")

        # Check solver status
        status = self.solver.get('status')
        if status != 0:
            pass
            #raise ValueError(f"HPIPM solver failed with status {status}")

        # Extract solution
        N_tot = sum(self.opts.horizon_lengths)
        x_trajectory = []
        u_trajectory = []
        for k in range(N_tot + 1):
            xk = self.qp_sol.get('x', k)
            x_trajectory.append(xk.copy())
            if k < N_tot:
                uk = self.qp_sol.get('u', k)
                u_trajectory.append(uk.copy())

        #stats = self.solver.get('time_ext')

        # Build final dictionary
        solution = {
            "status": "optimal" if status == 0 else "failed",
            "objective_value": 0.0,  # HPIPM doesn't directly return objective
            "x": x_trajectory,       # reduced-state trajectory
            "u": u_trajectory,
            "solve_time": solve_time
        }
        # Reconstruct the full-state trajectory
        solution["x_full"] = self.compute_full_trajectory(solution["x"])
        
        return solution

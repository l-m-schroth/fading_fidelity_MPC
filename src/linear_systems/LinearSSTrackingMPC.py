import numpy as np
import cvxpy as cp
import control as ct
from dataclasses import dataclass
#import matlab.engine
import os
from linear_systems.mor.balancedTruncationUnstable import balanced_truncation_unstable
import time

@dataclass
class LinearSSTrackingMPCOptions:
    """Options to configure the multi-fidelity MPC."""
    model_orders: list           # e.g., [n1, n2, ...]
    horizon_lengths: list        # e.g., [N1, N2, ...], sums to total horizon
    step_sizes: list             # e.g., [dt_0, dt_1, ..., dt_(N_tot-1)]
    Q: np.ndarray                # State tracking weight
    R: np.ndarray                # Control effort weight
    K: np.ndarray                # Terminal cost weight (symmetric PSD)
    ground_truth_ct_system: ct.StateSpace  # Continuous-time LTI system
    A_constr: np.ndarray
    b_constr: np.ndarray
    soft_constraints: bool = None
    lam: float = 0                # soft constraint cost scale

class LinearSSTrackingMPC:
    def __init__(self, opts: LinearSSTrackingMPCOptions):
        """
        Precompute everything needed for building and solving the MPC problem.
        """
        self.opts = opts
        
        # Basic checks
        assert len(self.opts.model_orders) == len(self.opts.horizon_lengths), (
            "The length of model_orders must match the length of horizon_lengths."
        )
        N_tot = sum(self.opts.horizon_lengths)
        assert len(self.opts.step_sizes) == N_tot, (
            "step_sizes must have a length equal to the total number of time steps."
        )

        # Parse dimensions from the ground truth system
        A_ct = self.opts.ground_truth_ct_system.A
        self.n_full = A_ct.shape[0]
        self.m = self.opts.ground_truth_ct_system.B.shape[1]

        # store matlab engine
        script_dir = os.path.dirname(os.path.abspath(__file__))  # Get script directory
        mor_path = os.path.join(script_dir, 'mor')  # Append 'mor' subdirectory
        #self.eng = matlab.engine.start_matlab()
        #self.eng.addpath(mor_path, nargout=0)

        # Precompute discrete-time reduced-order models:
        # We'll store them in a dictionary keyed by (n_red, dt_k).
        self.reduced_models = {}
        self._build_reduced_models()

        # Build optimization problem
        self.build_MPC_problem()
    
    def _build_reduced_models(self):
        """
        For each unique (model_order, dt) that appears in the problem,
        discretize the continuous-time system and then reduce to that order.
        """
        from math import isclose

        # Identify all pairs (model_order, dt) we need.
        # For chunk i in range(len(self.opts.model_orders)), we have N_i steps 
        # each with dt_k in step_sizes, so collect them all.
        needed_pairs = []
        dt_index = 0
        for i, n_red in enumerate(self.opts.model_orders):
            N_i = self.opts.horizon_lengths[i]
            for _ in range(N_i):
                # if dt_index >= len(self.opts.step_sizes)-1:
                #     break
                dt_k = self.opts.step_sizes[dt_index] # + 1 due to initial state
                needed_pairs.append((n_red, dt_k))
                dt_index += 1
        # add pair for initial state
        #needed_pairs.append((self.opts.model_orders[0], self.opts.step_sizes[0]))

        # Convert to a set to avoid re-computing duplicates
        needed_pairs = set(needed_pairs)

        # Build each discrete-time, reduced-order model
        if False:
            for (n_red, dt) in needed_pairs:
                sys_gt = ct.ss(self.opts.ground_truth_ct_system.A, self.opts.ground_truth_ct_system.B, np.eye(self.opts.ground_truth_ct_system.A.shape[0]), np.zeros((self.opts.ground_truth_ct_system.A.shape[0], self.opts.ground_truth_ct_system.B.shape[1])))
                sys_gt_discrete= ct.c2d(sys_gt, dt, method='zoh')


                # Convert Python arrays to MATLAB format
                Af_mat = matlab.double(sys_gt_discrete.A.tolist())
                Bf_mat = matlab.double(sys_gt_discrete.B.tolist())
                Cf_mat = matlab.double(np.eye(self.opts.ground_truth_ct_system.A.shape[0]).tolist())

                # Call MATLAB function
                A_red, B_red, C_red, W, V, S = self.eng.balancedTruncationUnstable(Af_mat, Bf_mat, Cf_mat, n_red, False, nargout=6)
                
                # Store in dictionary
                self.reduced_models[(n_red, dt)] = (np.array(A_red), np.array(B_red), np.array(C_red), np.array(W), np.array(V), np.array(S))
        else:
            for (n_red, dt) in needed_pairs:
                A_red_py, B_red_py, C_red_py, W_py, V_py, S_py = balanced_truncation_unstable(self.opts.ground_truth_ct_system.A, self.opts.ground_truth_ct_system.B, self.opts.ground_truth_ct_system.C, n_red, True)
                sys_reduced = ct.ss(A_red_py, B_red_py, C_red_py, np.zeros((C_red_py.shape[0], B_red_py.shape[1])))
                sys_reduced_discrete = ct.c2d(sys_reduced, dt, method='zoh')
                self.reduced_models[(n_red, dt)] = (sys_reduced_discrete.A, sys_reduced_discrete.B, sys_reduced_discrete.C, W_py, V_py, S_py)

    def build_MPC_problem(self):
        """
        Build the MPC problem
        """
        # ------------------ SETUP DIMENSIONS & VARIABLES ------------------
        N_tot = sum(self.opts.horizon_lengths)  # total steps

        # x[k] is the state at time step k (dimension depends on chunk_index(k))
        x_vars = []
        for k in range(N_tot):
            i = self.chunk_index(k) 
            n_red_k = self.opts.model_orders[i]
            x_vars.append(cp.Variable((n_red_k,), name=f"x_{k}"))
        x_vars.append(cp.Variable((self.opts.model_orders[-1],), name=f"x_{k}")) # terminal state

        # u[k] is the control at time step k
        u_vars = [cp.Variable((self.m,), name=f"u_{k}") for k in range(N_tot)]

        # Parameter for the initial condition in the *first* chunk.
        x0_param = cp.Parameter((self.opts.ground_truth_ct_system.A.shape[0],), value=None)

        # ------------------ BUILD COST & CONSTRAINTS ------------------
        constraints = []
        objective = 0.0

        # 2) For each time step, apply the appropriate discrete model and bridging
        for k in range(N_tot):
            i = self.chunk_index(k) # Identify which chunk k belongs to
            n_red_i = self.opts.model_orders[i]
            dt_k = self.opts.step_sizes[k]

            # Retrieve discrete-time system (A,B,C) for (n_red_i, dt_k)
            A_tilde, B_tilde, _, W, V, _ = self.reduced_models[(n_red_i, dt_k)]

            if k == 0:
                # 1) Initial condition
                constraints.append(x_vars[0] == W.T @ x0_param)
        
            # Build bridging matrix from chunk i dimension to chunk j dimension
            if k == N_tot - 1:
                Psi_ij = np.eye(self.opts.model_orders[-1]) # map to terminal state
            else:
                # Identify the dimension for x[k+1]
                j = self.chunk_index(k+1)
                n_red_j = self.opts.model_orders[j]
                dt_k_plus_1 = self.opts.step_sizes[k+1]
                # Retrieve next Model
                _, _, _, W_next, _, _ = self.reduced_models[(n_red_j, dt_k_plus_1)]
                Psi_ij = W_next.T @ V # matrix for brdging the model dimensions

            # Variables
            xk = x_vars[k]      # dimension n_red_i
            uk = u_vars[k]      # dimension m
            xk_plus = x_vars[k+1]  # dimension n_red_j

            # 2a) Dynamic constraint:
            #     x_{k+1} = Psi_ij * (A_i x_k + B_i u_k)
            constraints.append(
                xk_plus == Psi_ij @ A_tilde @ xk + Psi_ij @ B_tilde @ uk
            )

            # Add inequality constraints if A_constr and b_constr are not None
            if self.opts.A_constr is not None and self.opts.b_constr is not None:
                if self.opts.soft_constraints:
                    # Add convex soft constraint
                    # Define a slack variable (one per constraint)
                    s_u = cp.Variable(self.opts.b_constr.shape, nonneg=True, name=f"s_{k}")
                    s_l = cp.Variable(self.opts.b_constr.shape, nonneg=True, name=f"s_{k}")
                    # Replace the hard constraint with a softened version using the slack variable.
                    constraints.append(
                        self.opts.A_constr @ V @ xk <= self.opts.b_constr + s_u
                    )
                    constraints.append(
                        -self.opts.b_constr <= self.opts.A_constr @ V @ xk + s_l
                    )
                    # Add a penalty for using slack to the objective.
                    objective += dt_k * self.opts.lam * cp.sum(s_u)
                    objective += dt_k * self.opts.lam * cp.sum(s_l)
                else:
                    constraints.append(
                        self.opts.A_constr @ V @ xk <= self.opts.b_constr
                    )

            # 2b) Stage cost (scaled by dt_k):
            cost_stage = cp.quad_form(V @ xk, dt_k*self.opts.Q) + cp.quad_form(uk, dt_k*self.opts.R)
            objective += cost_stage

        # 3) Terminal cost on x[N_tot]
        dt_end = self.opts.step_sizes[-1]
        n_red_end = self.opts.model_orders[-1]
        _, _, _, _, V, _ = self.reduced_models[(n_red_end, dt_end)]
        x_terminal = x_vars[-1]
        objective += cp.quad_form(V @ x_terminal, self.opts.K) # weigh K with first dt, we typically select the K to be the discrete time (first dt) LQR 

        # ------------------ BUILD & SOLVE PROBLEM ------------------
        self.prob = cp.Problem(cp.Minimize(objective), constraints)
        self.x0_param = x0_param
        self.x_vars = x_vars
        self.u_vars = u_vars
    
    def solve(self, x0: np.ndarray):
        # Assign initial condition
        self.x0_param.value = x0
        
        start_time = time.time()
        self.prob.solve(eps_abs=1e-9, eps_rel=1e-9, max_iter = 10000)  # choose solver later
        solve_time = time.time() - start_time
        if self.prob.status not in [cp.OPTIMAL, cp.OPTIMAL_INACCURATE]:
            raise ValueError(f"Problem not solved to optimality. Status: {self.prob.status}")

        # ------------------ EXTRACT SOLUTION ------------------
        solution = {
            "status": self.prob.status,
            "objective_value": self.prob.value,
            "x": [x.value for x in self.x_vars],
            "u": [u.value for u in self.u_vars],
            "solve_time": solve_time
        }
        solution["x_full"] = self.compute_full_trajectory(solution["x"])

        return solution
    
    def chunk_index(self, k: int) -> int:
        """
        Returns the index of the chunk to which time step k belongs.
        For example, if horizon_lengths = [5, 5], steps 0..4 are chunk 0 and steps 5..9 are chunk 1.
        """
        s = 0
        for i, N_i in enumerate(self.opts.horizon_lengths):
            if k < s + N_i:
                return i
            s += N_i
        raise IndexError(f"Step k={k} is out of the valid range.")

    def compute_full_trajectory(self, reduced_trajectory: list) -> list:
        """
        Given a list of reduced state trajectories (one per time step), compute and return
        the corresponding full state trajectory for the original system. For a ROM,
        the full state is given by V @ x, where V is extracted from the corresponding
        reduced model used at that time step.
        """
        full_trajectory = []

        for k, x_red in enumerate(reduced_trajectory):
            if k == len(reduced_trajectory) - 1: # treat terminal state seperately
                dt_end = self.opts.step_sizes[-1]
                n_red_end = self.opts.model_orders[-1]
                _, _, _, _, V_end, _ = self.reduced_models[(n_red_end, dt_end)]
                full_trajectory.append(V_end @ x_red)
                break
            chunk = self.chunk_index(k) 
            step_size = self.opts.step_sizes[k]
            key = (self.opts.model_orders[chunk], step_size)
            # V is the 5th element in the tuple stored in self.reduced_models.
            V_mat = self.reduced_models[key][4]
            full_trajectory.append(V_mat @ x_red)
        return full_trajectory



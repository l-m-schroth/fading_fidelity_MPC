import numpy as np
import control as ct
from hpipm_python import *
from hpipm_python.common import *
from dataclasses import dataclass
from linear_systems.LinearSSTrackingMPC import LinearSSTrackingMPCOptions
import time

class LinearSSTrackingMPChpipmCondensing:
    """
    Single-model HPIPM-based MPC that directly discretizes the ground-truth
    system for each dt_k, then condenses out the states into a dense QP.
    Takes in LinearSSTrackingMPCOptions, where the discretization step sizes can be specified
    """
    def __init__(self, opts: LinearSSTrackingMPCOptions):
        self.opts = opts
        
        # Check single-chunk + full-order usage
        if len(opts.model_orders) != 1 or len(opts.horizon_lengths) != 1:
            raise ValueError(
                "Condensing requires exactly one chunk (single model_order, single horizon_length)."
            )
        self.N = opts.horizon_lengths[0]
        if len(opts.step_sizes) != self.N:
            raise ValueError("step_sizes must match horizon_length for single-chunk condensing.")
        
        self.n_full = opts.ground_truth_ct_system.A.shape[0]
        if opts.model_orders[0] != self.n_full:
            raise ValueError("Condensing must use the full-order model.")
        
        self.m = opts.ground_truth_ct_system.B.shape[1]

        # Discretize the ground-truth system for each dt_k
        self.A_list = []
        self.B_list = []
        sys_ct = self.opts.ground_truth_ct_system
        for dt in self.opts.step_sizes:
            sys_d = ct.c2d(sys_ct, dt, method='zoh')
            self.A_list.append(sys_d.A)
            self.B_list.append(sys_d.B)

        # Build condensed (A_bar, B_bar)
        self._build_condensed_dynamics()
        # Build the base dense QP (ignoring x0)
        self._build_dense_qp()

    def _build_condensed_dynamics(self):
        """Construct A_bar[k], B_bar[k] so x_{k+1} = A_bar[k]*x0 + B_bar[k]*z."""
        self.A_bar = []
        self.B_bar = []
        Phi_prev = np.eye(self.n_full)
        Psi_prev = np.zeros((self.n_full, self.m*self.N))

        for k in range(self.N):
            A_k = self.A_list[k]
            B_k = self.B_list[k]
            Phi_k = A_k @ Phi_prev
            Psi_k = A_k @ Psi_prev

            # Insert B_k into columns [k*m : (k+1)*m]
            blk = np.zeros_like(Psi_k)
            blk[:, k*self.m : (k+1)*self.m] = B_k
            Psi_k += blk

            self.A_bar.append(Phi_k)
            self.B_bar.append(Psi_k)

            Phi_prev = Phi_k
            Psi_prev = Psi_k

    def _build_dense_qp(self):
        """
        Create the 'base' dense QP structures that do not depend on x0.
        We'll incorporate x0 offsets in solve().
        """
        M = self.m * self.N  # dimension of z
        self.H_base = np.zeros((M, M))
        self.g_base = np.zeros(M)   # stays zero, ignoring x0 offset

        # Build cost: 0.5 * z^T H z + g^T z
        # Each stage cost => dt_k*( x_k^T Q x_k + u_k^T R u_k ), no extra factor
        for k in range(self.N):
            dt_k = self.opts.step_sizes[k]
            Qk = dt_k * self.opts.Q
            Rk = dt_k * self.opts.R

            if k == 0:
                x_k_matrix = np.zeros((self.n_full, M))
            else:
                x_k_matrix = self.B_bar[k-1]

            # Add Rk
            u_slice = slice(k*self.m, (k+1)*self.m)
            self.H_base[u_slice, u_slice] += Rk

            # Add Qk
            self.H_base += x_k_matrix.T @ Qk @ x_k_matrix

        # Terminal cost
        xN_matrix = self.B_bar[self.N - 1]
        self.H_base += xN_matrix.T @ self.opts.K @ xN_matrix

        # Now constraints ignoring x0
        # We want (for each stage) C z <= d
        C_list = []
        d_list = []
        if self.opts.A_constr is not None and self.opts.b_constr is not None:
            A_cons = self.opts.A_constr
            b_cons = 0.0*self.opts.b_constr # Set the constraint to zero, updated in solve
            for k in range(self.N):
                if k == 0:
                    x_k_matrix = np.zeros((self.n_full, M))
                else:
                    x_k_matrix = self.B_bar[k-1]
                C_k = A_cons @ x_k_matrix
                d_k = b_cons.copy()
                C_list.append(C_k)
                d_list.append(d_k)
            # terminal
            C_N = A_cons @ xN_matrix
            d_N = b_cons.copy()
            C_list.append(C_N)
            d_list.append(d_N)

        if len(C_list) > 0:
            self.C_base = np.vstack(C_list)
            self.d_base = np.hstack(d_list)
        else:
            self.C_base = np.zeros((0, M))
            self.d_base = np.zeros(0)

        # HPIPM dimension
        self.dim_dense = hpipm_dense_qp_dim()
        self.dim_dense.set('nv', M)
        self.dim_dense.set('ne', 0)
        self.dim_dense.set('nb', 0)
        self.dim_dense.set('ng', self.C_base.shape[0])

        # Build QP objects
        self.qp_dense = hpipm_dense_qp(self.dim_dense)
        self.qp_sol_dense = hpipm_dense_qp_sol(self.dim_dense)

        # Fill base H, C
        self.qp_dense.set('H', self.H_base)
        self.qp_dense.set('C', self.C_base)

        # We'll set 'g', 'ug', 'lg' at solve() time
        self.arg_dense = hpipm_dense_qp_solver_arg(self.dim_dense, 'robust')
        self.arg_dense.set('mu0', 1e4)
        self.arg_dense.set('iter_max', 200)
        self.arg_dense.set('tol_stat', 1e-6)
        self.arg_dense.set('tol_eq', 1e-6)
        self.arg_dense.set('tol_ineq', 1e-6)
        self.arg_dense.set('tol_comp', 1e-6)
        self.arg_dense.set('reg_prim', 1e-12)

        self.solver_dense = hpipm_dense_qp_solver(self.dim_dense, self.arg_dense)

    def solve(self, x0: np.ndarray):
        """
        Incorporate x0 in the cost offset g and constraints offset d, then solve.
        For now, only 'hard' constraints Cz <= d are handled. 
        """
        N = self.N
        M = self.m * N

        # 1) Build the linear offset g_offset from x0
        g_offset = np.zeros(M)
        for k in range(N):
            dt_k = self.opts.step_sizes[k]
            Qk = dt_k * self.opts.Q
            if k == 0:
                x_k_offset = x0
                x_k_matrix = np.zeros((self.n_full, M))
            else:
                x_k_offset = self.A_bar[k-1] @ x0
                x_k_matrix = self.B_bar[k-1]

            # derivative of 0.5 z^T H z => H z
            # offset => + x_k_matrix^T Qk x_k_offset
            g_offset += x_k_matrix.T @ Qk @ x_k_offset

        # Terminal offset
        xN_offset = self.A_bar[N-1] @ x0
        xN_matrix = self.B_bar[N-1]
        g_offset += xN_matrix.T @ self.opts.K @ xN_offset

        g_total = self.g_base + g_offset

        # 2) Build constraint offset => Cz <= d_total
        d_offset_list = []
        if self.opts.A_constr is not None and self.opts.b_constr is not None:
            A_con = self.opts.A_constr
            b_con = self.opts.b_constr
            for k in range(N):
                if k == 0:
                    x_k_off = x0
                else:
                    x_k_off = self.A_bar[k-1] @ x0
                d_k = b_con - A_con @ x_k_off
                d_offset_list.append(d_k)
            # terminal
            xN_off = self.A_bar[N-1] @ x0
            d_N = b_con - A_con @ xN_off
            d_offset_list.append(d_N)

        if len(d_offset_list) > 0:
            d_offset = np.concatenate(d_offset_list)
        else:
            d_offset = np.zeros(0)

        d_total = self.d_base + d_offset
        ng = self.C_base.shape[0]

        # 3) In HPIPM, to impose 'Cz <= d_total', we do:
        #    ug = d_total,  enable upper side with ug_mask=1
        #    lg = -∞ or disable lower side => we can set big negative or mask=0
        self.qp_dense.set('g', g_total)

        # Upper side
        self.qp_dense.set('ug', d_total)

        # If we want to enforce all these as "Cz <= d_total", we set mask=1 for each row
        ug_mask = np.ones(ng, dtype=int)
        self.qp_dense.set('ug_mask', ug_mask)

        # Lower side => disable or set to large negative
        # We'll do a big negative so that it's effectively -∞
        lg_val = -1e12 * np.ones(ng)
        self.qp_dense.set('lg', lg_val)
        # disable them => mask=0
        lg_mask = np.zeros(ng, dtype=int)
        self.qp_dense.set('lg_mask', lg_mask)

        # 4) Solve
        start_time = time.time()
        self.solver_dense.solve(self.qp_dense, self.qp_sol_dense)
        solve_time = time.time()- start_time

        status = self.solver_dense.get('status')
        if status != 0:
            print(f"[Condensing] HPIPM solver returned status={status}")

        # 5) Extract solution
        z_opt = self.qp_sol_dense.get('v')
        u_trajectory = []
        for k in range(N):
            u_k = z_opt[k*self.m : (k+1)*self.m]
            u_trajectory.append(u_k)

        return {
            "status": "optimal" if status == 0 else "failed",
            "objective_value": 0.0, 
            "u": u_trajectory,
            "solve_time": solve_time
        }


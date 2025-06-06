import numpy as np
import sys
import os
#import matlab.engine
import control as ct
from scipy.linalg import solve_discrete_are, solve_discrete_lyapunov
from scipy.optimize import fsolve
import numpy as np
from scipy.linalg import expm, solve
from scipy.integrate import fixed_quad

import numpy as np
from scipy.linalg import expm, solve
from scipy.integrate import fixed_quad, quad

import pickle
import os
import sys

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, parent_dir)  

# Get the parent directory
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
# mor_path = os.path.join(parent_dir, 'mor')  # Append 'mor' subdirectory

# Add parent directory to sys.path
sys.path.append(parent_dir)

from linear_systems.mor.balancedTruncationUnstable import balanced_truncation_unstable

def discrete_time_algebraic_ricatti_recursion(A_list, B_list, Q_list, R_list, P_terminal):
    """
    Computes the optimal controller and cost to go for the time varying LQR
    """
    K_list = []
    P_list = []
    P_t_plus_1 = P_terminal
    for A_t, B_t, Q_t, R_t in zip(list(reversed(A_list)), list(reversed(B_list)), list(reversed(Q_list)), list(reversed(R_list))):
        term = R_t + B_t.T @ P_t_plus_1 @ B_t
        K = -np.linalg.solve(term, B_t.T @ P_t_plus_1 @ A_t)
        P_t = Q_t + A_t.T @ P_t_plus_1 @ A_t - A_t.T @ P_t_plus_1 @ B_t @ np.linalg.inv(term) @ B_t.T @ P_t_plus_1 @ A_t
        P_t_plus_1 = P_t
        K_list.append(K)
        P_list.append(P_t)
    return list(reversed(K_list)), list(reversed(P_list))


def get_LQR_info(A_gt_cont, B_gt_cont, Q, R, model_orders, horizon_lengths, step_sizes, terminal_mode = 'full'):
    """
    Function for investigating LQR behaviour under model switching
    """

    A_list, B_list, Q_list, R_list = [], [], [], []
    reduced_models = build_reduced_models(A_gt_cont, B_gt_cont, model_orders, horizon_lengths, step_sizes)
    key_list = []

    # Initialize lists for inital state
    transition_idx = 0
    for i, n_red in enumerate(model_orders):
        N_i = horizon_lengths[i]
        for _ in range(N_i):
            dt = step_sizes[transition_idx]
            A, B, C, W, V, _ = reduced_models[(n_red, dt)]
            A_approx = V @ A @ W.T
            B_approx = V @ B
            A_list.append(A_approx)
            B_list.append(B_approx)
            Q_list.append(Q*dt)
            R_list.append(R*dt)
            key_list.append((n_red, dt))
            transition_idx += 1
    
    #print("A_list_lqr:", A_list)
    #print(key_list)
    # determine terminal costs
    if terminal_mode == 'full':
        sys_exact = ct.ss(A_gt_cont, B_gt_cont, np.eye(A_gt_cont.shape[0]), np.zeros((A_gt_cont.shape[0], B_gt_cont.shape[1])))
        dt = step_sizes[0]
        sys_exact_discrete = ct.c2d(sys_exact, dt, method='zoh')
        
        A_d = sys_exact_discrete.A
        B_d = sys_exact_discrete.B
        
        # Compute terminal cost using DARE
        P_terminal = solve_discrete_are(A_d, B_d, Q * dt, R * dt)

    elif terminal_mode == 'ROM':
        dt = step_sizes[0]
        A, B, _, W, V, _ = reduced_models[(model_orders[0], dt)]
        A_approx = V @ A @ W.T
        B_approx = V @ B
        
        # Compute terminal cost using DARE
        P_terminal = solve_discrete_are(A_approx, B_approx, Q * dt, R * dt)
    else:
        raise ValueError("Unknown terminal_mode. Choose 'full' or 'ROM'.")
    
    K_list, P_list = discrete_time_algebraic_ricatti_recursion(A_list, B_list, Q_list, R_list, P_terminal)

    def open_loop_plan(x0):
        ol_plan = [x0]
        x = x0
        for i in range(sum(horizon_lengths)):
            u = K_list[i] @ x
            x_next = A_list[i] @ x + B_list[i] @ u 
            x = x_next
            ol_plan.append(x)

        return ol_plan

    return K_list, P_list, open_loop_plan

def get_closed_loop_behaviour(x0, K, A, B, N):
    x_list = [x0]
    u_list = []
    x = x0
    for i in range(N):
        u = K @ x
        x_next = A @ x + B @ u 
        x = x_next
        x_list.append(x)
        u_list.append(u)
    return x_list, u_list

def compute_P_Lyapunov(A, B, K, Q, R):
    A_cl = A + B @ K
    Q_cl = Q + K.T @ R @ K
    
    # Solve the discrete-time Lyapunov equation:
    # A_cl.T @ P @ A_cl - P + Q_cl = 0
    P = solve_discrete_lyapunov(A_cl.T, Q_cl)
    return P

def build_reduced_models(A_gt_cont, B_gt_cont, model_orders, horizon_lengths, step_sizes):
    """
    For each unique (model_order, dt) that appears in the problem,
    discretize the continuous-time system and then reduce to that order.
    """

    # Identify all pairs (model_order, dt) we need.
    # For chunk i in range(len(self.opts.model_orders)), we have N_i steps 
    # each with dt_k in step_sizes, so collect them all.
    needed_pairs = []
    dt_index = 0
    for i, n_red in enumerate(model_orders):
        N_i = horizon_lengths[i]
        for _ in range(N_i):
            dt_k = step_sizes[dt_index] # + 1 due to initial state
            needed_pairs.append((n_red, dt_k))
            dt_index += 1

    # add pair for initial state and ground truth discretized system
    needed_pairs.append((model_orders[0], step_sizes[0]))

    # Convert to a set to avoid re-computing duplicates
    needed_pairs = set(needed_pairs)

    # Build each discrete-time, reduced-order model
    reduced_models = {}
    for (n_red, dt) in needed_pairs:

        sys_gt = ct.ss(A_gt_cont, B_gt_cont, np.eye(A_gt_cont.shape[0]), np.zeros((A_gt_cont.shape[0], B_gt_cont.shape[1])))
        sys_gt_discrete= ct.c2d(sys_gt, dt, method='zoh')


        # #Convert Python arrays to MATLAB format
        # Af_mat = matlab.double(sys_gt_discrete.A.tolist())
        # Bf_mat = matlab.double(sys_gt_discrete.B.tolist())
        # Cf_mat = matlab.double(np.eye(A_gt_cont.shape[0]).tolist())

        # #Call MATLAB function
        # A_red, B_red, C_red, W, V, S = eng.balancedTruncationUnstable(Af_mat, Bf_mat, Cf_mat, n_red, False, nargout=6)
        # Call Python implementation of balanced truncation
        A_red, B_red, C_red, W, V, S = balanced_truncation_unstable(sys_gt_discrete.A, sys_gt_discrete.B, sys_gt_discrete.C, n_red, False)

        # Store in dictionary
        reduced_models[(n_red, dt)] = (np.array(A_red), np.array(B_red), np.array(C_red), np.array(W), np.array(V), np.array(S))

    return reduced_models

def get_number_optimization_vars(model_schedule, horizon_lengths, u_dim, condesing):
    N_vars = 0 
    if not condesing:
        # State variables
        N_vars += model_schedule[0] # initial state
        for i, horizon_lens in enumerate(horizon_lengths):
            N_vars += model_schedule[i]*horizon_lens
    # Input variables
    N_vars += u_dim*(sum(horizon_lengths))
    return N_vars

def get_exponential_rates(N_steps_dt_increase_min, N, step=2):
    rates = []
    steps = []
    n_steps = N 
    while n_steps >= N_steps_dt_increase_min:
        def eq(alpha):
            return (1 - alpha**n_steps)/(1 - alpha) - N
        alpha_sol = fsolve(eq, 1.1)[0]
        rates.append(alpha_sol)
        steps.append(n_steps)
        n_steps -= step
    return rates, [[step] for step in steps]

def get_model_schedules(max_order_init, model_orders_fixed_model, absolute_min_model_order, N_models_max, step):
    max_order = max_order_init
    model_schedules = []
    while max_order >= max_order_init//2:
        diff = max_order - absolute_min_model_order
        for scale in [0.1, 0.25, 0.5, 0.75, 1]:
            min_order = int(max_order - scale*diff)
            for n_models in range(2,N_models_max+1):
                if scale == 0.1:
                    n_models_x = 2
                else:
                    n_models_x = n_models
                model_orders = []
                for numerator in range(n_models_x):
                    model_orders.append(int(np.ceil((max_order - numerator/(n_models-1)*(max_order - min_order)))))
                model_schedules.append(model_orders)               
        max_order -= step
    model_schedules.extend([[model_order] for model_order in model_orders_fixed_model])
    return model_schedules

def divide_evenly(N, num_chunks):
    quotient = N // num_chunks
    remainder = N % num_chunks
    return [quotient + 1 if i < remainder else quotient for i in range(num_chunks)]

def compute_Q_tilde_1_matrices(A, B, Q, dt, n_quad=10):
    """
    Computes the three matrices (Q1_1, Q1_2, Q1_3) that make up Q_tilde_1 under
    zero-order hold, following the same 'element-by-element' integration style
    from your old code.
    
    The formulas (from your LaTeX) are:
        Q1_1 = ∫[0..dt] e^(A^T t) * Q * e^(A t) dt
        Q1_2 = B^T ∫[0..dt] M(t)^T * Q * M(t) dt * B
        Q1_3 = 2 * ∫[0..dt] e^(A^T t) * Q * M(t) dt * B
        
    where M(t) = ∫[0..t] e^(A (t - τ)) dτ = A^(-1)( e^(A t) - I )  (if A is invertible).
    
    Parameters
    ----------
    A : np.ndarray, shape (n, n)
        Continuous-time system matrix.
    B : np.ndarray, shape (n, m)
        Input matrix.
    Q : np.ndarray, shape (n, n)
        State weighting matrix.
    dt : float
        Zero-order hold sampling period (time step).
    n_quad : int, optional
        Number of quadrature points for `fixed_quad`.
        
    Returns
    -------
    Q1_1 : (n, n) ndarray
    Q1_2 : (n, n) ndarray
    Q1_3 : (n, n) ndarray
    """
    n = A.shape[0]

    # -- Optional check if A is invertible (same style as your old code) --
    detA = np.linalg.det(A)
    if abs(detA) < 1e-15:
        raise ValueError(
            "Matrix A is (close to) singular. "
            "Consider a different approach or using np.linalg.pinv(A)."
        )

    # Helper: matrix exponential
    def eA(t):
        return expm(A * t)

    # Helper: M(t) = ∫[0..t] e^(A (t-τ)) dτ = A^(-1)(e^(A t) - I)
    def M_single_t(t):
        return solve(A, eA(t) - np.eye(n))

    # ----------------------------------------------------------------------
    #  Q1_1 = ∫[0..dt] e^(A^T t) * Q * e^(A t) dt
    # ----------------------------------------------------------------------
    def Q1_1_single_t(t):
        eAt = eA(t)
        return eAt.T @ Q @ eAt

    def integrand_element_Q1_1(t_array, i, j):
        """Return the (i,j) entry of Q1_1_single_t(t) for each t in t_array."""
        if np.isscalar(t_array):
            return Q1_1_single_t(t_array)[i, j]
        else:
            # If t_array is a vector, compute each then collect into a 1D array
            results = []
            for t in t_array:
                val_ij = Q1_1_single_t(t)[i, j]
                results.append(val_ij)
            return np.array(results)

    Q1_1 = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            Q1_1[i, j], _ = fixed_quad(
                integrand_element_Q1_1,  # integrand
                0, dt,                  # integration limits
                args=(i, j),
                n=n_quad
            )

    # ----------------------------------------------------------------------
    #  Q1_2 = B^T ( ∫[0..dt] M(t)^T Q M(t) dt ) B
    #        We first integrate the inner matrix over t, then multiply by B^T, B.
    # ----------------------------------------------------------------------
    def Q1_2_single_t(t):
        Mt = M_single_t(t)
        return Mt.T @ Q @ Mt

    def integrand_element_Q1_2(t_array, i, j):
        """Return the (i,j) entry of M(t)^T Q M(t)."""
        if np.isscalar(t_array):
            return Q1_2_single_t(t_array)[i, j]
        else:
            results = []
            for t in t_array:
                val_ij = Q1_2_single_t(t)[i, j]
                results.append(val_ij)
            return np.array(results)

    Q1_2_inner = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            val, _ = fixed_quad(
                integrand_element_Q1_2,
                0, dt,
                args=(i, j),
                n=n_quad
            )
            Q1_2_inner[i, j] = val

    Q1_2 = B.T @ Q1_2_inner @ B

    # ----------------------------------------------------------------------
    #  Q1_3 = 2 * ( ∫[0..dt] e^(A^T t) Q M(t) dt ) B
    # ----------------------------------------------------------------------
    def Q1_3_single_t(t):
        eAt = eA(t)
        Mt = M_single_t(t)
        return eAt.T @ Q @ Mt  # shape: (n, n)

    def integrand_element_Q1_3(t_array, i, j):
        """Return the (i,j) entry of e^(A^T t) Q M(t)."""
        if np.isscalar(t_array):
            return Q1_3_single_t(t_array)[i, j]
        else:
            results = []
            for t in t_array:
                val_ij = Q1_3_single_t(t)[i, j]
                results.append(val_ij)
            return np.array(results)

    Q1_3_inner = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            val, _ = fixed_quad(
                integrand_element_Q1_3,
                0, dt,
                args=(i, j),
                n=n_quad
            )
            Q1_3_inner[i, j] = val

    Q1_3 = 2.0 * Q1_3_inner @ B

    return Q1_1, Q1_2, Q1_3


def compute_Q_tilde(K, Q1_1, Q1_2, Q1_3, R, T):
    Q_tilde_1 = Q1_1 + K.T @ Q1_2 @ K + Q1_3 @ K
    Q_tilde_2 = T * K.T @ R @ K
    return Q_tilde_1 + Q_tilde_2

import numpy as np
from scipy.linalg import expm, solve
from scipy.integrate import fixed_quad

def compute_Q_tilde_1_matrices_efficient(A, B, Q, dt, n_quad=10):

    # ----------------------------------------------------------------------
    # Preliminary checks
    # ----------------------------------------------------------------------
    n = A.shape[0]
    m = B.shape[1]

    detA = np.linalg.det(A)
    if abs(detA) < 1e-15:
        raise ValueError("Matrix A is singular or nearly singular. "
                         "Consider using a pseudo-inverse or numerical integration of M(t).")

    # ----------------------------------------------------------------------
    # Caching: store e^(A t) and M(t) for each t in the quadrature nodes
    # ----------------------------------------------------------------------
    eA_cache = {}
    M_cache = {}

    def eA(t):
        """Compute/cached e^(A t)."""
        if t not in eA_cache:
            eA_cache[t] = expm(A * t)
        return eA_cache[t]

    def M(t):
        """
        M(t) = ∫[0..t] e^(A (t - tau)) d tau = A^(-1)( e^(A t) - I ).
        Assumes A is invertible.
        """
        if t not in M_cache:
            M_cache[t] = solve(A, eA(t) - np.eye(n))
        return M_cache[t]

    # ----------------------------------------------------------------------
    # Q1_1 = ∫[0..dt] e^(A^T t) Q e^(A t) dt
    # ----------------------------------------------------------------------
    Q1_1 = np.zeros((n, n))

    def integrand_Q1_1(t):
        EA = eA(t)
        return EA.T @ Q @ EA  # shape (n,n)

    def integrand_element_Q1_1(t_array, i, j):
        """Returns the (i,j)-element of e^(A^T t) Q e^(A t) for each t in t_array."""
        # If a single scalar t, just compute it once
        if np.isscalar(t_array):
            return integrand_Q1_1(t_array)[i, j]
        else:
            # Otherwise, for an array of times, build a result array
            results = []
            for t in t_array:
                val_ij = integrand_Q1_1(t)[i, j]
                results.append(val_ij)
            return np.array(results)

    for i in range(n):
        for j in range(n):
            val, _ = quad(integrand_element_Q1_1, 0, dt, args=(i, j))
            Q1_1[i, j] = val

    # ----------------------------------------------------------------------
    # Q1_2 = B^T ( ∫[0..dt] M(t)^T Q M(t) dt ) B
    #   We'll build the integral piece ∫[0..dt] [M(t)^T Q M(t)] dt first
    # ----------------------------------------------------------------------
    Q1_2_inner = np.zeros((n, n))

    def integrand_Q1_2(t):
        Mt = M(t)
        return Mt.T @ Q @ Mt

    def integrand_element_Q1_2(t_array, i, j):
        if np.isscalar(t_array):
            return integrand_Q1_2(t_array)[i, j]
        else:
            results = []
            for t in t_array:
                val_ij = integrand_Q1_2(t)[i, j]
                results.append(val_ij)
            return np.array(results)

    for i in range(n):
        for j in range(n):
            val, _ = quad(integrand_element_Q1_2, 0, dt, args=(i, j))
            Q1_2_inner[i, j] = val

    Q1_2 = B.T @ Q1_2_inner @ B

    # ----------------------------------------------------------------------
    # Q1_3 = 2 * ∫[0..dt] e^(A^T t) Q M(t) dt  B
    # ----------------------------------------------------------------------
    Q1_3_inner = np.zeros((n, n))

    def integrand_Q1_3(t):
        EA = eA(t)
        Mt = M(t)
        return EA.T @ Q @ Mt

    def integrand_element_Q1_3(t_array, i, j):
        if np.isscalar(t_array):
            return integrand_Q1_3(t_array)[i, j]
        else:
            results = []
            for t in t_array:
                val_ij = integrand_Q1_3(t)[i, j]
                results.append(val_ij)
            return np.array(results)

    for i in range(n):
        for j in range(n):
            val, _ = quad(integrand_element_Q1_3, 0, dt, args=(i, j))
            Q1_3_inner[i, j] = val

    Q1_3 = 2.0 * Q1_3_inner @ B

    return Q1_1, Q1_2, Q1_3

def get_random_linear_system(n, m, l, C_eye=True, discrete=True):
    # Generate stable and controllable system matrices
    while True:
        print("Creating random linear system")

        # Generate a stable, discrete-time, random state-space system
        if discrete:
            sys = ct.drss(n, l, m, dt=True)
        else:
            sys = ct.rss(n, l, m)

        A = sys.A
        B = sys.B
        if C_eye and n==l:
            C = np.eye(n)
        else:
            C = sys.C

        eigs = np.linalg.eigvals(A)
        if discrete:
            print("largest norm of eigenvalue of A:", np.max(np.abs(eigs)))
        else:
            print("largest norm of eigenvalue of A:", np.max(np.real(eigs)))

        # Check if (A, B) is controllable
        ctrb_matrix = ct.ctrb(A, B)
        obsv_matrix = ct.obsv(A, C)
        controllability = np.linalg.matrix_rank(ctrb_matrix, tol=1e-5) == n
        observability = np.linalg.matrix_rank(obsv_matrix, tol=1e-5) == n
        if controllability and observability:
            break

    return A, B, C

def generate_augmented_linear_system(n, n_unstable, m, l, C_eye=True, discrete=True):
    "generate a random linear system with stable and unstable modes"
    A_stable, B_stable, C_stable = get_random_linear_system(n - n_unstable, m, l, C_eye=C_eye, discrete=discrete)
    
    # Create unstable modes (e.g., eigenvalues uniformly chosen between 0.1 and 1.0)
    if not discrete: # continuous case_
        unstable_modes = np.diag(np.random.uniform(0.1, 1.0, size=n_unstable))
    else:
        # For discrete systems, choose eigenvalues greater than 1 (e.g., in [1.1, 2.0])
        unstable_modes = np.diag(np.random.uniform(1.1, 2.0, size=n_unstable))
    
    # Augment A: block-diagonal combination of the stable and unstable parts
    A_full = np.block([
        [A_stable,              np.zeros((n - n_unstable, n_unstable))],
        [np.zeros((n_unstable, n - n_unstable)), unstable_modes]
    ])
    while True:
        # Extend B: stack the stable B with a randomly generated unstable B part.
        B_unstable = 0.1*np.random.randn(n_unstable, m)
        B_full = np.vstack([B_stable, B_unstable])
        
        # Extend C:
        if C_eye:
            # Force C to be the identity matrix of size n.
            # (For this to make sense, you typically need l == n.)
            C_full = np.eye(n)
        else:
            # Extend the stable C (size: (l, n - n_unstable)) with a random block (size: (l, n_unstable))
            C_unstable = 0.1*np.random.randn(l, n_unstable)
            C_full = np.hstack([C_stable, C_unstable])
        
        # Check controllability and observability of the augmented system
        ctrb_matrix = ct.ctrb(A_full, B_full)
        obsv_matrix = ct.obsv(A_full, C_full)
        
        controllability = np.linalg.matrix_rank(ctrb_matrix, tol=1e-5) == n
        observability = np.linalg.matrix_rank(obsv_matrix, tol=1e-5) == n
        if controllability and observability:
            # Found a system that is both controllable and observable.
            break
        else:
            print("Generated system not controllable/observable. Regenerating...")
    
    # Define D as zero (assuming no direct feedthrough term)
    D_full = np.zeros((C_full.shape[0], B_full.shape[1]))
    return A_full, B_full, C_full, D_full

def closed_loop_simulation_2(simulation_frequency, control_frequency, mpc, x0, duration, error_scale=0):
    # Determine the simulation and control time steps
    dt_sim = 1.0 / simulation_frequency
    dt_ctrl = 1.0 / control_frequency

    # Discretize the continuous ground truth system at the simulation rate
    sys_d = ct.c2d(mpc.opts.ground_truth_ct_system, dt_sim, method='zoh')
    A_d = np.array(sys_d.A)
    B_d = np.array(sys_d.B)
    C_d = np.array(sys_d.C)
    
    # Extract cost matrices
    Q = mpc.opts.Q
    R = mpc.opts.R

    # Initialize trajectories and metrics
    states = [x0]
    outputs = [C_d @ x0]
    controls = []          # Store control updates (only when computed)
    costs = []             # Store stage cost at control update instants
    costs_slack = []       # cost of the soft constraint
    costs_total = []
    comp_times = []        # Store computation times at control update instants
    
    x_current = x0.copy()

    # Compute total simulation steps and steps between control updates
    total_steps = int(np.ceil(duration / dt_sim))
    steps_per_control = int(np.ceil(dt_ctrl / dt_sim))
    
    # Initialize the control input (could be zeros or an initial guess)
    u_current = np.zeros(B_d.shape[1])
    
    for i in range(total_steps):
        # Only update the control at intervals of dt_ctrl (zero order hold in between)
        if i % steps_per_control == 0:
            # Solve the MPC problem from the current state
            result = mpc.solve(x_current)
            comp_times.append(result["solve_time"])
            
            # Extract the first control input from the MPC solution
            u_seq = result["u"]
            u_current = np.array(u_seq[0]).flatten()  # ensure a flat vector
            controls.append(u_current)
            
        # Compute and record stage cost at the control update
        stage_cost = x_current.T @ Q @ x_current + u_current.T @ R @ u_current 
        # Add soft constraint penalty if applicable
        if mpc.opts.soft_constraints:
            s = np.maximum(0, mpc.opts.A_constr @ x_current - mpc.opts.b_constr)
            stage_cost_slack = mpc.opts.lam * np.sum(s)
        else: 
            stage_cost_slack = 0
        costs_slack.append(stage_cost_slack)
        costs.append(stage_cost)
        costs_total.append(stage_cost + stage_cost_slack)
        
        # Propagate the state using the discretized model and current control input
        x_next = A_d @ x_current + B_d @ u_current
        # Add process noise scaled by dt_sim if error_scale > 0
        x_current = x_next + error_scale * dt_sim * np.random.randn(*x_next.shape)
        
        states.append(x_current)
        outputs.append(C_d @ x_current)
    
    closed_loop_traj = {
        "states": states,
        "outputs": outputs,
        "controls": controls,
        "costs": costs,
        "costs_slack":costs_slack,
        "costs_total":costs_total,
        "computation_times": comp_times
    }
    
    return closed_loop_traj

def generate_exponential_step_sizes(initial_step: float, rate: float, num_steps: int, plot: bool = False) -> list:
    """
    Generate a list of exponentially increasing step sizes.

    Parameters:
        initial_step (float): The initial step size.
        rate (float): The exponential growth rate (e.g., 1.05 for a 5% increase per step).
        num_steps (int): The total number of step sizes to generate.
        plot (bool): If True, display a plot of the step size development.

    Returns:
        list: A list of exponentially increasing step sizes.
    """
    step_sizes = [initial_step * (rate ** i) for i in range(num_steps)]

    if plot:
        import matplotlib.pyplot as plt
        plt.figure()
        plt.plot(range(num_steps), step_sizes, marker='o')
        plt.title("Exponential Step Sizes")
        plt.xlabel("Step Index")
        plt.ylabel("Step Size")
        plt.grid(True)
        plt.show()

    return step_sizes


import numpy as np
from scipy.linalg import solve_continuous_lyapunov, solve_discrete_lyapunov, cholesky

def balanced_truncation(Af, Bf, Cf, k, continuous=True):
    """
    Perform balanced truncation on a stable continuous (or discrete) LTI system.
    
    The system is assumed to be controllable (Af, Bf) and observable (Af, Cf).
    This function implements the balancing-free square-root algorithm.
    
    Inputs:
        Af, Bf, Cf : np.ndarray
            Full-order system matrices.
        k : int
            Reduced model order.
        continuous : bool, optional
            True for continuous-time systems, False for discrete-time.
    
    Returns:
        A_red, B_red, C_red : np.ndarray
            Reduced-order system matrices.
        W, V : np.ndarray
            Petrov-Galerkin projection matrices (satisfying W'*V = I).
        s_vals : np.ndarray
            Hankel singular values (vector of singular values).
    """
    # Compute the controllability Gramian factor.
    eps = 1e-6 # need to add small epsilon to cholesky, otherwise it tends to fail for matrices that are close to singular
    if continuous:
        # Solve A*X + X*A' + B*B' = 0 for X, then compute its Cholesky factor.
        X = solve_continuous_lyapunov(Af, -Bf @ Bf.T)
        # cholesky returns a lower triangular factor L such that L @ L.T = X.
        # In MATLAB the function "lyapchol" returns a lower triangular factor which is then transposed.
        U = cholesky(X + eps*np.eye(X.shape[0]), lower=True)
        
        # Compute the observability Gramian factor.
        Y = solve_continuous_lyapunov(Af.T, -Cf.T @ Cf)
        L = cholesky(Y + eps*np.eye(Y.shape[0]), lower=True)
    else:
        # For discrete-time: X = A*X*A' + B*B'
        X = solve_discrete_lyapunov(Af, Bf @ Bf.T)
        U = cholesky(X + eps*np.eye(X.shape[0]), lower=True)
        
        Y = solve_discrete_lyapunov(Af.T, Cf.T @ Cf)
        L = cholesky(Y + eps*np.eye(Y.shape[0]), lower=True)

    # Perform SVD 
    W_svd, s_vals, Vh = np.linalg.svd(U.T @ L, full_matrices=False)
    V_svd = Vh.T

    # Perform QR factorizations:
    # X is the orthogonal factor from U*W_svd
    X_qr, _ = np.linalg.qr(U @ W_svd)
    # Y is the orthogonal factor from L_mat*V_svd
    Y_qr, _ = np.linalg.qr(L @ V_svd)

    # Truncate to the first k columns
    Xk = X_qr[:, :k]
    Yk = Y_qr[:, :k]

    # Compute projection matrices.
    # In MATLAB: W = Yk*inv(Xk'*Yk); V = Xk;
    W_proj = Yk @ np.linalg.inv(Xk.T @ Yk)
    V_proj = Xk

    # Compute reduced system matrices.
    A_red = W_proj.T @ Af @ V_proj
    B_red = W_proj.T @ Bf
    C_red = Cf @ V_proj

    return A_red, B_red, C_red, W_proj, V_proj, s_vals
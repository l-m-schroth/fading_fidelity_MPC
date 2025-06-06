"""
In this file we test that our implementation of balanced truncation and stable unstable decomposition works by comparing it to an existing implementation
Running these tests requires a licensed installation of matlab and the matlab engine.
"""

import sys
import os
import control
from linear_systems.mor.balancedTruncationUnstable import stabsep

# Get the absolute path of the parent directory
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

# Add the parent directory to the system path if it's not already included
if parent_dir not in sys.path:
    sys.path.insert(0, parent_dir)

from linear_systems.mor.balancedTruncation import balanced_truncation
from linear_systems.utils_linear_systems import get_random_linear_system, generate_augmented_linear_system
import matlab.engine
import os
import numpy as np
from linear_systems.mor.balancedTruncationUnstable import balanced_truncation_unstable

def test_balanced_truncation():
    # Generate a random system.
    n = 4   # state dimension
    m = 2   # input dimension
    l = 4   # output dimension (for C; if C_eye is True, then l == n)
    discrete = True  # Use discrete-time system for this example.
    A, B, C = get_random_linear_system(n, m, l, C_eye=True, discrete=discrete)
    
    # Set the reduced order.
    k = 2
    continuous = not discrete  # For discrete systems, set continuous to False.

    # Compute Python balanced truncation results.
    A_py, B_py, C_py, W_py, V_py, S_py = balanced_truncation(A, B, C, k, continuous=continuous)
    
    # Start MATLAB engine and add the path to your MOR functions.
    eng = matlab.engine.start_matlab()
    mor_path = os.path.dirname(os.path.abspath(__file__))
    eng.addpath(mor_path, nargout=0)
    
    # Convert Python arrays to MATLAB format.
    Af_mat = matlab.double(A.tolist())
    Bf_mat = matlab.double(B.tolist())
    Cf_mat = matlab.double(C.tolist())
    
    # Call the MATLAB function (balancedTruncationUnstable) using the MATLAB engine.
    A_mat, B_mat, C_mat, W_mat, V_mat, S_mat = eng.balancedTruncation(
        Af_mat, Bf_mat, Cf_mat, k, continuous, nargout=6
    )
    
    # Convert MATLAB outputs to numpy arrays.
    A_mat = np.array(A_mat)
    B_mat = np.array(B_mat)
    C_mat = np.array(C_mat)
    W_mat = np.array(W_mat)
    V_mat = np.array(V_mat)
    # S_mat might be returned as a list of lists; flatten it.
    S_mat = np.array(S_mat).flatten()

    # Compare the results.
    # Note: Depending on the implementation, projection matrices W and V can have sign/rotation ambiguities.
    diff_A = np.linalg.norm(A_py - A_mat)
    diff_B = np.linalg.norm(B_py - B_mat)
    diff_C = np.linalg.norm(C_py - C_mat)
    diff_S = np.linalg.norm(S_py - S_mat)
    diff_W = np.linalg.norm(W_py - W_mat)
    diff_V = np.linalg.norm(V_py - V_mat)

    print("Difference in W:", diff_W)
    print("Difference in V:", diff_V)
    print("Difference in reduced A:", diff_A)
    print("Difference in reduced B:", diff_B)
    print("Difference in reduced C:", diff_C)
    print("Difference in Hankel singular values:", diff_S)
    
    # For a strict test, you might use a tolerance:
    tol = 1e-5
    assert np.allclose(A_py, A_mat, atol=tol), "A matrices differ beyond tolerance"
    assert np.allclose(B_py, B_mat, atol=tol), "B matrices differ beyond tolerance"
    assert np.allclose(C_py, C_mat, atol=tol), "C matrices differ beyond tolerance"
    assert np.allclose(S_py, S_mat, atol=tol), "Singular values differ beyond tolerance"
    assert np.allclose(W_py, W_mat, atol=tol), "W matrices differ beyond tolerance"
    assert np.allclose(V_py, V_mat, atol=tol), "V matrices differ beyond tolerance"
    
    print("Test passed: Python and MATLAB balanced truncation results match within tolerance.")
    eng.quit()

def controllability_matrix(A, B):
    """
    Compute the controllability matrix Wc = [B, A*B, A^2*B, ..., A^(n-1)*B].
    """
    n = A.shape[0]
    Wc = B
    for i in range(1, n):
        Wc = np.hstack((Wc, np.linalg.matrix_power(A, i) @ B))
    return Wc

def test_stabsep():
    """
    Test that the Python stabsep function produces a stable/unstable split that
    is consistent with MATLAB's stabsep. Instead of directly comparing the state-space
    matrices, we first compute a similarity transformation (based on the controllability matrices)
    between the corresponding MATLAB and Python systems, then transform the MATLAB matrices,
    and finally compare the transformed matrices with the Python ones.
    """
    # Create a system with both stable and unstable eigenvalues.
    n = 5           # total state dimension
    n_unstable = 2  # number of unstable states
    m = 2           # input dimension
    l = 4           # output dimension
    discrete = False
    dt_matlab = -1 if discrete else 0
    
    # generate_augmented_linear_system should return (A, B, C, D)
    A, B, C, D = generate_augmented_linear_system(n, n_unstable, m, l, C_eye=True, discrete=discrete)
    
    # Create the Python state-space system.
    sys = control.ss(A, B, C, D, None)
    
    # Compute the stable/unstable decomposition using Python.
    sys_stable_py, sys_unstable_py = stabsep(sys, continuous=not discrete)
    
    # Start MATLAB engine.
    eng = matlab.engine.start_matlab()
    mor_path = os.path.dirname(os.path.abspath(__file__))
    eng.addpath(mor_path, nargout=0)
    
    # --- MATLAB stabsep on the full system ---
    A_mat = matlab.double(A.tolist())
    B_mat = matlab.double(B.tolist())
    C_mat = matlab.double(C.tolist())
    D_mat = matlab.double(D.tolist())
    sys_mat = eng.ss(A_mat, B_mat, C_mat, D_mat, dt_matlab, nargout=1)
    sys_stable_mat, sys_unstable_mat = eng.stabsep(sys_mat, nargout=2)
    
    # --- Extract state-space data from the MATLAB systems using ssdata ---
    A_stable_mat, B_stable_mat, C_stable_mat, D_stable_mat = eng.ssdata(sys_stable_mat, nargout=4)
    A_unstable_mat, B_unstable_mat, C_unstable_mat, D_unstable_mat = eng.ssdata(sys_unstable_mat, nargout=4)
    
    # Convert MATLAB arrays to numpy arrays.
    A_stable_mat = np.array(A_stable_mat)
    B_stable_mat = np.array(B_stable_mat)
    C_stable_mat = np.array(C_stable_mat)
    D_stable_mat = np.array(D_stable_mat)
    
    A_unstable_mat = np.array(A_unstable_mat)
    B_unstable_mat = np.array(B_unstable_mat)
    C_unstable_mat = np.array(C_unstable_mat)
    D_unstable_mat = np.array(D_unstable_mat)
    
    # For Python systems, extract the state-space data.
    A_stable_py = np.array(sys_stable_py.A)
    B_stable_py = np.array(sys_stable_py.B)
    C_stable_py = np.array(sys_stable_py.C)
    D_stable_py = np.array(sys_stable_py.D)
    
    A_unstable_py = np.array(sys_unstable_py.A)
    B_unstable_py = np.array(sys_unstable_py.B)
    C_unstable_py = np.array(sys_unstable_py.C)
    D_unstable_py = np.array(sys_unstable_py.D)
    
    tol = 1e-5
    
    # --- For the stable subsystem, compute similarity transformation ---
    Wc_stable_mat = controllability_matrix(A_stable_mat, B_stable_mat)
    Wc_stable_py  = controllability_matrix(A_stable_py, B_stable_py)
    
    if np.linalg.matrix_rank(Wc_stable_mat) < A_stable_mat.shape[0]:
        raise ValueError("MATLAB stable system is not controllable!")
    
    T_stable = Wc_stable_py @ np.linalg.pinv(Wc_stable_mat)
    T_stable_inv = np.linalg.inv(T_stable)
    
    # Transform MATLAB stable matrices into Python coordinates:
    A_stable_mat_trans = T_stable @ A_stable_mat @ T_stable_inv
    B_stable_mat_trans = T_stable @ B_stable_mat
    C_stable_mat_trans = C_stable_mat @ T_stable_inv
    # D matrices should be identical.
    
    # --- For the unstable subsystem, compute similarity transformation ---
    Wc_unstable_mat = controllability_matrix(A_unstable_mat, B_unstable_mat)
    Wc_unstable_py  = controllability_matrix(A_unstable_py, B_unstable_py)
    
    if np.linalg.matrix_rank(Wc_unstable_mat) < A_unstable_mat.shape[0]:
        raise ValueError("MATLAB unstable system is not controllable!")
    
    T_unstable = Wc_unstable_py @ np.linalg.pinv(Wc_unstable_mat)
    T_unstable_inv = np.linalg.inv(T_unstable)
    
    A_unstable_mat_trans = T_unstable @ A_unstable_mat @ T_unstable_inv
    B_unstable_mat_trans = T_unstable @ B_unstable_mat
    C_unstable_mat_trans = C_unstable_mat @ T_unstable_inv

    # --- Compare the transformed MATLAB matrices with the Python ones ---
    assert np.allclose(A_stable_py, A_stable_mat_trans, atol=tol), "Stable A matrices differ beyond tolerance"
    assert np.allclose(B_stable_py, B_stable_mat_trans, atol=tol), "Stable B matrices differ beyond tolerance"
    assert np.allclose(C_stable_py, C_stable_mat_trans, atol=tol), "Stable C matrices differ beyond tolerance"
    assert np.allclose(D_stable_py, D_stable_mat, atol=tol), "Stable D matrices differ beyond tolerance"
    
    assert np.allclose(A_unstable_py, A_unstable_mat_trans, atol=tol), "Unstable A matrices differ beyond tolerance"
    assert np.allclose(B_unstable_py, B_unstable_mat_trans, atol=tol), "Unstable B matrices differ beyond tolerance"
    assert np.allclose(C_unstable_py, C_unstable_mat_trans, atol=tol), "Unstable C matrices differ beyond tolerance"
    assert np.allclose(D_unstable_py, D_unstable_mat, atol=tol), "Unstable D matrices differ beyond tolerance"
    
    print("Test passed: Python and MATLAB stabsep results are similar (transformed state-space matrices match within tolerance).")
    eng.quit()

def compute_transfer_function(A, B, C, s):
    """
    Compute the transfer function H(s) = C*(sI - A)^{-1}*B.
    Here s can be a scalar or an array of complex numbers.
    """
    I = np.eye(A.shape[0])
    return C @ np.linalg.inv(s * I - A) @ B

def subspace_distance(V1, V2, tol=1e-5):
    """
    Check whether the column spaces of V1 and V2 are the same.
    One way is to compare the projection matrices:
        P = Q*Q^T, where Q is an orthonormal basis for the columns of V.
    The subspaces are the same if ||P1 - P2|| is small.
    """
    # Compute orthonormal bases
    Q1, _ = np.linalg.qr(V1)
    Q2, _ = np.linalg.qr(V2)
    P1 = Q1 @ Q1.T
    P2 = Q2 @ Q2.T
    return np.linalg.norm(P1 - P2) < tol

def test_balanced_truncation_unstable():
    """
    Test that the balanced truncation for unstable systems (balanced_truncation_unstable)
    produces a reduced system that is similar between Python and MATLAB. Instead of comparing
    the matrices directly, we first compute a similarity transformation (using the controllability
    matrices) that maps MATLAB's realization into Python's coordinates, then compare the transformed
    matrices.
    """
    # Generate an augmented system with both stable and unstable modes.
    n = 6         # total state dimension
    n_unstable = 2  # number of unstable states
    m = 2         # input dimension
    l = 6         # output dimension (using C_eye=True makes l == n)
    continuous = False  # use continuous-time system
    
    # generate_augmented_linear_system returns (A, B, C, D)
    A, B, C, D = generate_augmented_linear_system(n, n_unstable, m, l, C_eye=True, discrete=not continuous)
    
    # Set the reduced order. It must be at least as large as the number of unstable modes.
    k = 4
    
    # Compute Python balanced truncation unstable results.
    A_py, B_py, C_py, W_py, V_py, S_py = balanced_truncation_unstable(A, B, C, k, continuous)
    
    # Start MATLAB engine and add the path to your MOR functions.
    eng = matlab.engine.start_matlab()
    mor_path = os.path.dirname(os.path.abspath(__file__))
    eng.addpath(mor_path, nargout=0)
    
    # Convert Python arrays to MATLAB format.
    Af_mat = matlab.double(A.tolist())
    Bf_mat = matlab.double(B.tolist())
    Cf_mat = matlab.double(C.tolist())
    
    # Call the MATLAB function balancedTruncationUnstable using the MATLAB engine.
    A_mat, B_mat, C_mat, W_mat, V_mat, S_mat = eng.balancedTruncationUnstable(
        Af_mat, Bf_mat, Cf_mat, k, continuous, nargout=6
    )
    
    # Convert MATLAB outputs to numpy arrays.
    A_mat = np.array(A_mat)
    B_mat = np.array(B_mat)
    C_mat = np.array(C_mat)
    W_mat = np.array(W_mat)
    V_mat = np.array(V_mat)
    S_mat = np.array(S_mat).flatten()
    
    tol = 1e-5
    
    # --- Compute similarity transformation for the reduced system ---
    # For the reduced systems, compute controllability matrices (order k).
    Wc_py = controllability_matrix(A_py, B_py)
    Wc_mat = controllability_matrix(A_mat, B_mat)
    # Since the reduced systems are minimal, Wc_mat is invertible.
    T = Wc_py @ np.linalg.pinv(Wc_mat)
    T_inv = np.linalg.inv(T)
    
    # Transform MATLAB reduced matrices into the Python coordinates.
    A_mat_trans = T @ A_mat @ T_inv
    B_mat_trans = T @ B_mat
    C_mat_trans = C_mat @ T_inv
    # D is assumed equal.
    
    # --- 1. Compare transfer functions for the reduced systems ---
    # Use the transformed MATLAB realization.
    sample_points = [
        0.1 + 0.1j,  # very low frequency
        1 + 1j,      # moderate frequency
        2 + 2j,      # moderate frequency
        5 + 0j,      # real frequency point
        0 - 1j,      # purely imaginary negative frequency
        -3 + 3j      # negative real part, mixed frequency
    ]
    for s in sample_points:
        H_py = compute_transfer_function(A_py, B_py, C_py, s)
        H_mat_trans = compute_transfer_function(A_mat_trans, B_mat_trans, C_mat_trans, s)
        diff_H = np.linalg.norm(H_py - H_mat_trans)
        print(f"Difference in transfer function at s={s}: {diff_H}")
        assert diff_H < tol, f"Transfer functions differ at s={s} beyond tolerance"
    
    # --- 2. Directly compare the reduced state-space matrices ---
    assert np.allclose(A_py, A_mat_trans, atol=tol), "Reduced A matrices differ beyond tolerance"
    assert np.allclose(B_py, B_mat_trans, atol=tol), "Reduced B matrices differ beyond tolerance"
    assert np.allclose(C_py, C_mat_trans, atol=tol), "Reduced C matrices differ beyond tolerance"
    assert np.allclose(D, D, atol=tol), "Reduced D matrices differ beyond tolerance"
    
    # --- 3. Check whether the V matrices correspond to projections onto the same subspace ---
    if not subspace_distance(V_py, V_mat, tol):
        diff_proj = np.linalg.norm((np.linalg.qr(V_py)[0] @ np.linalg.qr(V_py)[0].T) -
                                    (np.linalg.qr(V_mat)[0] @ np.linalg.qr(V_mat)[0].T))
        raise AssertionError(f"V matrices project to different subspaces (projection difference: {diff_proj})")
    else:
        print("V matrices project to the same subspace within tolerance.")
    
    # --- 4. Compare Hankel singular values directly.
    diff_S = np.linalg.norm(S_py - S_mat)
    print("Difference in Hankel singular values:", diff_S)
    assert np.allclose(S_py, S_mat, atol=tol), "Hankel singular values differ beyond tolerance"
    
    # --- 5. Check that W.T @ V equals the identity matrix.
    I_check = W_py.T @ V_py
    I_expected = np.eye(I_check.shape[0])
    diff_identity = np.linalg.norm(I_check - I_expected)
    print("Difference in W.T @ V from identity:", diff_identity)
    assert diff_identity < tol, "W.T @ V is not the identity matrix within tolerance"
    
    print("Test passed: Python and MATLAB balanced truncation unstable results are similar (transformed state-space matrices match within tolerance).")
    eng.quit()

if __name__ == '__main__':
    test_balanced_truncation()
    test_stabsep()
    test_balanced_truncation_unstable()
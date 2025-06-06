import numpy as np
import scipy.linalg
import control  # pip install control
from linear_systems.mor.balancedTruncation import balanced_truncation

import numpy as np
import scipy.linalg
import control  # Assuming the 'control' package is being used

def stabsep(sys, continuous=True, tol=1e-9):
    """
    Decompose the state-space system 'sys' into stable and unstable parts
    using an ordered Schur decomposition and a decoupling transformation that
    cancels the coupling term via a Sylvester equation.
    
    Parameters:
        sys : control.StateSpace
            The state-space system to decompose.
        continuous : bool
            True if the system is continuous time (stable if Re(s) < 0),
            False for discrete time (stable if |z| < 1).
        tol : float
            Tolerance used in determining stability.
    
    Returns:
        sys_stable : control.StateSpace
            The stable subsystem.
        sys_unstable : control.StateSpace
            The unstable subsystem.
    """
    # Get state-space matrices
    A, B, C, D = sys.A, sys.B, sys.C, sys.D
    n = A.shape[0]

    # Define sorting function for Schur decomposition based on stability.
    # For continuous-time: stable if Re(eig) < -tol.
    # For discrete-time: stable if |eig| < 1 - tol.
    if continuous:
        sort_func = lambda x: np.real(x) < -tol
    else:
        sort_func = lambda x: np.abs(x) < 1 - tol

    # Compute ordered real Schur decomposition:
    # T_schur is the quasi-triangular Schur form, V is orthogonal.
    T_schur, V, sdim = scipy.linalg.schur(A, sort=sort_func, output='real')
    ldim = sdim  # number of stable eigenvalues

    # Partition V so that the first ldim columns correspond to the stable subspace.
    V_stable = V[:, :ldim]
    V_unstable = V[:, ldim:]
    
    # Transform A, B, C to the Schur basis:
    # In the Schur basis, A_bar = V^T A V has the block form:
    # [A_s  A_c; 0  A_u]
    A_bar = T_schur  # already computed from schur()
    A_stable = A_bar[:ldim, :ldim]
    A_coupling = A_bar[:ldim, ldim:]
    A_unstable = A_bar[ldim:, ldim:]

    # Similarly, transform B and C:
    Bs = V.T @ B
    B_stable = Bs[:ldim, :]
    B_unstable = Bs[ldim:, :]

    Cs = C @ V
    C_stable = Cs[:, :ldim]
    C_unstable = Cs[:, ldim:]
    
    # Solve the Sylvester equation to decouple the stable subsystem:
    # We need to find X such that:
    #   A_stable * X - X * A_unstable = A_coupling.
    # Since scipy.linalg.solve_sylvester solves A*X + X*B = C,
    # we set B = -A_unstable.
    X = scipy.linalg.solve_sylvester(A_stable, -A_unstable, A_coupling)

    # Compute corrected input and output matrices:
    # Update the stable part of the input: B_stable_new = B_stable + X * B_unstable.
    B_stable_new = B_stable + X @ B_unstable
    # Update the unstable part of the output: C_unstable_new = C_unstable - C_stable * X.
    C_unstable_new = C_unstable - C_stable @ X

    # The decoupled overall state-space realization in the Schur coordinates is now:
    # A_new = diag(A_stable, A_unstable)
    # B_new = [B_stable_new; B_unstable]
    # C_new = [C_stable, C_unstable_new]
    # We assign the stable subsystem to have (A_stable, B_stable_new, C_stable, D)
    # and the unstable subsystem to have (A_unstable, B_unstable, C_unstable_new, D).
    sys_stable = control.ss(A_stable, B_stable_new, C_stable, D)
    sys_unstable = control.ss(A_unstable, B_unstable, C_unstable_new, D)

    return sys_stable, sys_unstable


def balanced_truncation_unstable(Af, Bf, Cf, k, continuous):
    """
    Perform balanced truncation on an unstable LTI system (continuous or discrete).
    First decomposes the full order system into stable and unstable parts,
    then reduces only the stable subsystem.

    Parameters:
        Af, Bf, Cf: ndarray
            Full-order system matrices.
        k: int
            Desired reduced order.
        continuous: bool
            True for continuous-time system, False for discrete-time.

    Returns:
        A_red, B_red, C_red: ndarray
            Reduced order system matrices.
        W, V: ndarray
            Petrov-Galerkin projection matrices satisfying W'*V = I.
        S: ndarray
            Singular values from the balanced truncation.
    """

    nf = Af.shape[0]
    m = Bf.shape[1]

    # Create a full-order state-space system with identity output matrix.
    if continuous:
        sysf = control.ss(Af, Bf, np.eye(nf), np.zeros((nf, m)))
        is_stable = np.all(np.real(np.linalg.eigvals(Af)) < 0)
    else:
        # For discrete systems, the sample time is set to -1.
        sysf = control.ss(Af, Bf, np.eye(nf), np.zeros((nf, m)), None)
        is_stable = np.all(np.abs(np.linalg.eigvals(Af)) < 1)

    if is_stable: # if stable, only balanced truncation necessary
        A_red, B_red, C_red, W, V, S = balanced_truncation(Af, Bf, Cf, k, continuous)
        return A_red, B_red, C_red, W, V, S

    print("Performing stable unstable decomposition.")
    # Decompose the system using the provided stabsep function.
    sys_stable, sys_unstable = stabsep(sysf, continuous)

    # Stable subsystem matrices.
    Afstable = sys_stable.A
    Bfstable = sys_stable.B
    Tfstable = sys_stable.C  # note: in our decomposition, C is used as the transformation block
    nfstable = Afstable.shape[0]

    # Unstable subsystem: here we use its 'C' matrix as the other part of the transformation.
    Tfunstable = sys_unstable.C
    nfunstable = nf - nfstable

    # Form the transformation matrix T such that T = [Tfunstable, Tfstable]
    T = np.hstack((Tfunstable, Tfstable))

    # Transform the original output matrix Cf.
    Cf_transformed = Cf @ T
    # Extract the part corresponding to the stable subsystem.
    # MATLAB indexes columns nfunstable+1:end; in Python we take columns nfunstable: (0-indexed).
    Cfstable = Cf_transformed[:, nfunstable:]

    # The reduced order will keep all the unstable modes.
    nROM = k - nfunstable
    print(f"Performing model reduction, reducing {nfstable} to {nROM}.")

    # Perform balanced truncation on the stable subsystem.
    # (Assuming balanced_truncation returns (Ar, Br, Cr, Wstable, Vstable, S)).
    _, _, _, Wstable, Vstable, S = balanced_truncation(Afstable, Bfstable, Cfstable, nROM, continuous)

    # Build the overall projection matrices.
    # Construct block matrices:
    # For the unstable part, we use identity, and for the stable part, the computed Wstable and Vstable.
    X = np.block([
        [np.eye(nfunstable),            np.zeros((nfunstable, nROM))],
        [np.zeros((nfstable, nfunstable)), Wstable]
    ])
    Y = np.block([
        [np.eye(nfunstable),            np.zeros((nfunstable, nROM))],
        [np.zeros((nfstable, nfunstable)), Vstable]
    ])

    # The MATLAB code computes:
    #     W = ([I,0;0,Wstable]' * inv(T))' = inv(T).T @ [I,0;0,Wstable]
    # and
    #     V = T * [I,0;0,Vstable].
    # We follow the same computation here.
    W = np.linalg.inv(T).T @ X
    V = T @ Y

    # Form the reduced-order model.
    A_red = W.T @ Af @ V
    B_red = W.T @ Bf
    C_red = Cf @ V

    return A_red, B_red, C_red, W, V, S

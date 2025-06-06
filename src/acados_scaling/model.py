from acados_template import AcadosModel
from casadi import SX, vertcat
import casadi as ca
import numpy as np
from typing import Tuple
import control

# define seed for reproducibility
np.random.seed(0)

def generate_random_linear_system(n: int, m: int) -> Tuple[ca.SX, ca.SX]:
    np.random.seed(1) # fix random seed for reproducibility
    system = control.rss(n, n, m, dt=True)
    A, B = system.A, system.B
    return ca.SX(A), ca.SX(B)

def export_model(n: int, m: int) -> AcadosModel:
    # generate random linear system
    Ak, Bk = generate_random_linear_system(n, m)

    # set up states & controls
    x = SX.sym('x', n)
    u = SX.sym('u', m)

    # xdot
    xdot = SX.sym('xdot', n)

    model = AcadosModel()

    model.disc_dyn_expr = Ak @ x + Bk @ u

    model.x = x
    model.xdot = xdot
    model.u = u
    model.name = 'random_linear_dynamics'
    
    return model


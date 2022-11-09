import numpy as np
from scipy.linalg import eigvalsh
import scipy.sparse as spsp
import scipy.sparse.linalg as spspla
from gputils import incomplete_cholesky

def condition_number(A):
    l = abs(eigvalsh(A))
    return max(l) / min(l)


def conjugate_gradient(A, b, x0=None, tol=1e-7):

    # Get the initial approximate solution
    if x0 is None:
        x = np.zeros(b.size)
    else:
        x = x0

    r = b - A @ x
    p = r
    k = 0

    while np.linalg.norm(r) > tol:
        alpha = (r @ r) / (p @ A @ p)

        x = x + alpha * p

        r_new = r - alpha * A @ p

        if np.linalg.norm(r_new) <= tol:
            break

        beta = (r_new @ r_new) / (r @ r)

        p = r_new + beta * p

        r = r_new

        k += 1

    return x


def preconditioned_conjugate_gradient(A, b, x0=None, tol=1e-7, P=None, L=None, get_history=False):

    history = []

    # Get the initial approximate solution
    if x0 is None:
        x = np.zeros(b.size)
    else:
        x = x0

    if get_history:
        history.append(x)

    # Get the preconditioner matrix if given
    if P is None:
        pass
    elif P == 'diag':
        P_inv = spsp.diags(1/A.diagonal(), format='csr')
        L = None
    elif P == 'ichol':
        P_inv = None
        if L is None:
            L = incomplete_cholesky(A)
    else:
        raise ValueError("P has to be either None, 'ichol' or 'diag'")

    def get_z(r):
        if P is None:
            return r
        elif P_inv is not None:
            return P_inv @ r
        elif L is not None:
            Lr = spspla.spsolve(L, r)
            LTLr = spspla.spsolve(L.T, Lr)
            return LTLr
        else:
            assert False, 'Either P is None or P_inv is calculated, or L is calculated'

    r = b - A @ x
    z = get_z(r)
    p = z
    k = 0

    while np.linalg.norm(r) > tol:
        alpha = (r @ z) / (p @ A @ p)

        x = x + alpha * p

        if get_history:
            history.append(x)

        r_new = r - alpha * A @ p

        if np.linalg.norm(r_new) <= tol:
            break

        z_new = get_z(r_new)

        beta = (r_new @ z_new) / (r @ z)

        p = z_new + beta * p

        r = r_new
        z = z_new

        k += 1

    if get_history:
        return history
    else:
        return x

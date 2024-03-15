import numpy as np
import scipy.sparse as spsp
from sksparse import cholmod as cm


def get_reorder(A, ordering_method=None):
    if ordering_method is None:
        ordering_method = "amd"

    # Transpose to convert the matrix from csr to csc
    chol = cm.analyze(A.T, ordering_method=ordering_method)
    reorder = chol.P()

    N = len(reorder)
    P = spsp.csr_array((np.ones(N), reorder, np.arange(N + 1)), shape=(N, N))

    return P


def reorder_vector(b, P):
    return P @ b


def rev_reorder_vector(b, P):
    return P.T @ b


def reorder_matrix(A, P):
    return P @ A @ P.T


def rev_reorder_matrix(A, P):
    return P.T @ A @ P

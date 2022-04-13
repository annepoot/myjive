"""
This module is intended to easily add functions to the GaussianModule / GaussianModel
It should be deprecated in the end, but helps set things up more easily.
"""

import numpy as np
import warnings

def get_phis(N_obs, N_mesh):

    # Create an empty array
    phi = np.zeros((N_mesh, N_obs))
    phi_sub = np.zeros((N_mesh, N_obs))

    # Get the indices of the nodes that will be observed
    vec_obs = np.rint(np.linspace(0, N_mesh-1, N_obs))
    vec_obs = vec_obs.astype(int)

    for i in range(N_obs):

        phi[vec_obs[i], i] = 1
        phi_sub[vec_obs[i], i] = 1

    j = 0
    j_old = 0

    for i in range(N_obs-1):
        while j < N_mesh:
            if np.isclose(phi[j,i+1], 1):
                phi[j_old:j+1,i] = np.linspace(1, 0, j-j_old+1)
                phi[j_old:j+1,i+1] = np.linspace(0, 1, j-j_old+1)
                j_old = j
                break

            j = j+1

    return phi, phi_sub

# Define the exponential covariance function
def cov_exp(x1, x2, **kwargs):

    # Get the relevant kwargs
    l = kwargs.get("l", 1)
    sigma_f = kwargs.get("sigma_f", 1)

    # Compute the covariance
    cov = sigma_f**2 * np.exp(-abs(x1-x2)/l)

    return cov

# Define the squared exponential covariance function
def cov_sq_exp(x1, x2, **kwargs):

    # Get the relevant kwargs
    l = kwargs.get("l", 1)
    sigma_f = kwargs.get("sigma_f", 1)

    # Compute the covariance
    cov = sigma_f**2 * np.exp(-0.5 * ((x1-x2)/l)**2)

    return cov

def cov_sq_exp_semicorr(x1, x2, **kwargs):

    cov = np.zeros(x1.shape)
    mid = kwargs.get("L", 1) / 2

    for i in range(x1.shape[0]):
        for j in range(x1.shape[1]):
            x1_item, x2_item = x1[i,j], x2[i,j]

            if x1_item < mid and x2_item < mid:
                cov[i,j] = cov_sq_exp(x1_item, x2_item, **kwargs)
            elif x1_item > mid and x2_item > mid:
                cov[i,j] = cov_sq_exp(x1_item, x2_item, **kwargs)
            elif np.isclose(x1_item, mid) and np.isclose(x2_item, mid):
                cov[i,j] = cov_sq_exp(x1_item, x2_item, **kwargs)

    return cov

def cov_poly_block(x1, x2, **kwargs):

    cov = np.zeros(x1.shape)
    mid = 0
    C = 1
    D = 1

    for i in range(x1.shape[0]):
        for j in range(x1.shape[1]):
            x1_item, x2_item = x1[i,j], x2[i,j]

            if x1_item < mid and x2_item < mid:
                cov[i,j] = 0.25 * x1_item**2 * x2_item**2 + C * x1_item * x2_item + D
            elif x1_item > mid and x2_item > mid:
                cov[i,j] = 0.25 * x1_item**2 * x2_item**2 + C * x1_item * x2_item + D
            elif np.isclose(x1_item, x2_item):
                cov[i,j] = 0.25 * x1_item**2 * x2_item**2 + C * x1_item * x2_item + D
            else:
                cov[i,j] = C * x1_item * x2_item + D

    return cov

def cov_block(x1, x2, **kwargs):

    cov = np.zeros(x1.shape)
    L = kwargs.get("L", 1)

    for i in range(x1.shape[0]):
        for j in range(x1.shape[1]):
            x1_item, x2_item = x1[i,j], x2[i,j]

            if np.isclose(x1_item, 0):
                if np.isclose(x2_item, 0):
                    cov[i,j] = 1
            elif np.isclose(x1_item, L):
                if np.isclose(x2_item, L):
                    cov[i,j] = 1
            else:
                if not np.isclose(x2_item, 0) and not np.isclose(x2_item, L):
                    cov[i,j] = 1

    return cov

def kernel(X1, X2, **kwargs):

    # Get the kwargs
    X2_mesh, X1_mesh = np.meshgrid(X2, X1)

    K = kwargs['cov_func'](X1_mesh, X2_mesh, **kwargs)

    # Check if we need to ensure positive-definite ness
    if kwargs.get('ensure_PD', False):

        # Get the smallest eigenvalue, and check if it is negative and close to 0
        eig = min(np.real(np.linalg.eigvalsh(K)))

        if np.isclose(eig,0) and eig < 0:

            # If so, add a tiny bit of noise to the covariance matrix
            noise = 10 * abs(eig)
            warnings.warn('Adding i.i.d. noise with magnitude {:.4e} to ensure positive definiteness of the covariance matrix'.format(noise))
            K += noise * np.identity(K.shape[0])

            # Get the smallest eigenvalue, and check if it is negative and close to 0
            eig = min(np.real(np.linalg.eig(K)[0]))

    return K


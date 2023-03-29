import sys
sys.path.append('../../')

import numpy as np
import scipy.sparse.linalg as spspla
from time import perf_counter
import jive.util.proputils as pu
from jive.app import main
from jive.solver.constrainer import Constrainer

import solverutils as su

import jive.solver.jit.cholesky as nbchol
import jive.solver.util.reorder as reord

props = pu.parse_file('beam.pro')

densities = ['coarse', 'medium', 'fine', 'fine2', 'fine3', 'fine4']

for density in densities:
    props['init']['mesh']['file'] = 'meshes/beam_' + density + '.msh'
    globdat = main.jive(props)

    K = globdat['matrix0']
    f = globdat['extForce']
    u = globdat['state0']
    c = globdat['constraints']

    conman = Constrainer(c, K)
    Kc = conman.get_output_matrix()
    fc = conman.get_rhs(f)

    N = K.shape[0]

    Kc_full = Kc.toarray()

    P = reord.get_reorder(Kc)

    p2 = su.get_reorder(Kc)
    P2 = su.get_reorder_matrix(p2)
    assert np.allclose(P.toarray(), P2.toarray())

    print(N)

    print('Performing reorder')
    time = -perf_counter()
    Kr = reord.reorder_matrix(Kc, P)
    time += perf_counter()
    print(time, '\n')

    print('Performing rev_reorder')
    time = -perf_counter()
    Krr = reord.rev_reorder_matrix(Kc, P)
    time += perf_counter()
    print(time, '\n')

    fr = reord.reorder_vector(fc, P)

    uc = spspla.spsolve(Kc, fc)
    assert np.allclose(u, uc)

    ur = spspla.spsolve(Kr, fr)
    assert np.allclose(u, reord.rev_reorder_vector(ur, P))

    if N < 2000:
        print('Performing np.linalg.cholesky')
        time = -perf_counter()
        L0 = np.linalg.cholesky(Kc_full)
        time += perf_counter()
        print(time, '\n')

        print('Performing nbchol.sparse_cholesky')
        time = -perf_counter()
        L1 = nbchol.sparse_cholesky(Kc)
        time += perf_counter()
        print(time, '\n')

    print('Performing nbchol.sparse_cholesky (reordered)')
    time = -perf_counter()
    L2 = nbchol.sparse_cholesky(Kr)
    time += perf_counter()
    print(time, '\n')

    print('Performing nbchol.incomplete_cholesky')
    time = -perf_counter()
    L3 = nbchol.incomplete_cholesky(Kc)
    time += perf_counter()
    print(time, '\n')

    print('Performing nbchol.incomplete_cholesky (reordered)')
    time = -perf_counter()
    L4 = nbchol.incomplete_cholesky(Kr)
    time += perf_counter()
    print(time, '\n')

    if N < 2000:
        assert np.allclose(L0, L1.toarray())
        assert np.allclose((L1 @ L1.T).toarray(), Kc.toarray())
    assert np.allclose((L2 @ L2.T).toarray(), Kr.toarray())

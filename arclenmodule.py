import warnings
import numpy as np
import scipy.sparse as sparse
import scipy.sparse.linalg as linalg

from numpy.linalg import norm as norm

from names import GlobNames as gn
from names import ParamNames as pn
from names import Actions as act

from module import *
from constrainer import Constrainer

NSTEPS = 'nsteps'
ITERMAX = 'itermax'
TOLERANCE = 'tolerance'
BETA = 'beta'
DL = 'dl'


class ArclenModule(Module):
    def init(self, props, globdat):
        myprops = props[self._name]
        self._step = 0
        self._nsteps = int(myprops.get(NSTEPS, 1))
        self._itermax = int(myprops.get(ITERMAX, 100))
        self._tolerance = float(myprops.get(TOLERANCE, 1e-6))
        self._beta = float(myprops.get(BETA, 0.1))
        self._dl = float(myprops.get(DL, 0.02))

        model = globdat[gn.MODEL]
        dc = globdat[gn.DOFSPACE].dof_count()
        self._fext0 = np.zeros(dc)
        self._fhat = np.zeros(dc)
        params = {pn.EXTFORCE: self._fext0, pn.UNITFORCE: self._fhat}
        model.take_action(act.GETEXTFORCE, params, globdat)
        model.take_action(act.GETUNITFORCE, params, globdat)

    def run(self, globdat):
        dc = globdat[gn.DOFSPACE].dof_count()
        model = globdat[gn.MODEL]

        globdat[gn.TIMESTEP] = self._step
        print('Running time step', self._step)

        if self._step == 0:
            globdat[gn.STATE0] = np.zeros(dc)
            globdat[gn.LAMBDA] = 0.
            duOld = np.zeros(dc)
        else:
            duOld = globdat[gn.STATE0] - globdat[gn.OLDSTATE0]

        globdat[gn.OLDSTATE0] = np.copy(globdat[gn.STATE0])
        K = np.zeros((dc, dc))
        fint = np.zeros(dc)
        fhat = self._fhat
        c = Constrainer()

        params = {
            pn.MATRIX0: K,
            pn.INTFORCE: fint,
            pn.CONSTRAINTS: c,
        }

        # Initialize first iteration
        iteration = 0

        # Advance to next time step
        model.take_action(act.ADVANCE, params, globdat)

        # Assemble K
        model.take_action(act.GETMATRIX0, params, globdat)

        # Assemble fext0
        fext0 = fhat * globdat[gn.LAMBDA] + self._fext0

        # Get constraints
        model.take_action(act.GETCONSTRAINTS, params, globdat)
        cdofs, cvals = c.get_constraints()
        fdofs = [i for i in range(dc) if i not in cdofs]

        # Solve system
        duII = solveSys(K, self._fhat, c)

        if self._step < 1:
            lsign = 1
        else:  # Choose solution with smallest angle wrt previous solution
            lsign = np.sign(np.dot(duOld, fhat) * np.dot(duII, fhat))

        beta2F = self._beta ** 2 * np.dot(fhat, fhat)
        Dlam = lsign * self._dl / np.sqrt(np.dot(duII, duII) + beta2F)
        lbeta2F = Dlam * beta2F

        Du0 = Dlam * duII
        globdat[gn.STATE0] += Du0

        # Reference values to check convergence
        ref = max(Dlam * np.linalg.norm(fext0), 1)
        rel = 1.

        # Initialize iteration loop
        while rel > self._tolerance and iteration < self._itermax:
            iteration += 1
            params[pn.MATRIX0] = np.zeros((dc, dc))
            model.take_action(act.GETMATRIX0, params, globdat)
            params[pn.INTFORCE] = np.zeros(dc)
            model.take_action(act.GETINTFORCE, params, globdat)
            fext = fext0 + Dlam * fhat
            r = fext - params[pn.INTFORCE]

            # Solve arclength system of equations(sherman-morrison method):
            du1 = solveSys(params[pn.MATRIX0], r, c)
            du2 = solveSys(params[pn.MATRIX0], -fhat, c)
            vdI = np.dot(Du0, du1)
            vdIIk = np.dot(Du0, du2 - lbeta2F)
            du = du1 - vdI / vdIIk * du2
            dl = vdI / vdIIk

            globdat[gn.STATE0] += du
            Dlam += dl
            rel = np.linalg.norm(r[np.ix_(fdofs)]) / ref
            print('Iteration %i, relative residual norm: %.4e' % (iteration, rel))

        # Alert if not convergence
        if rel > self._tolerance:
            if rel > 1:
                raise RuntimeError('Divergence in time step %i' % self._step)
            else:
                warnings.warn('No convergence in time step %i' % self._step)

        globdat[gn.LAMBDA] += Dlam

        self._step += 1

        if self._step >= self._nsteps:
            return 'exit'
        else:
            return 'ok'

    def shutdown(self, globdat):
        pass


def solveSys(K, f, c):
    Kc, fc = c.constrain(K, f)
    smat = sparse.csr_matrix(Kc)
    u = linalg.spsolve(smat, fc)
    return u


def declare(factory):
    factory.declare_module('Arclen', ArclenModule)

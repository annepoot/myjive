import warnings
import numpy as np
import scipy.sparse as sparse
import scipy.sparse.linalg as linalg

from jive.names import GlobNames as gn
from jive.names import ParamNames as pn
from jive.names import Actions as act

from jive.app import Module
from jive.solver import Constrainer

NSTEPS = "nsteps"
ITERMAX = "itermax"
TOLERANCE = "tolerance"
BETA = "beta"
DL = "dl"

__all__ = ["ArclenModule"]


class ArclenModule(Module):
    def init(self, props, globdat):
        myprops = props[self._name]
        self._step = 0
        self._nsteps = int(myprops.get(NSTEPS, 1))
        self._itermax = int(myprops.get(ITERMAX, 100))
        self._tolerance = float(myprops.get(TOLERANCE, 1e-6))
        self._beta = float(myprops.get(BETA, 0.1))
        self._dl = float(myprops.get(DL, 0.02))
        globdat[gn.ACCEPTED] = True

        self._models = globdat[gn.MODELS]
        dc = globdat[gn.DOFSPACE].dof_count()
        self._fext0 = np.zeros(dc)
        self._fhat = np.zeros(dc)
        params = {pn.EXTFORCE: self._fext0, pn.UNITFORCE: self._fhat}

        for model in self.get_relevant_models("GETEXTFORCE", self._models):
            model.GETEXTFORCE(params, globdat)

        for model in self.get_relevant_models("GETUNITFORCE", self._models):
            model.GETUNITFORCE(params, globdat)

    def run(self, globdat):
        dc = globdat[gn.DOFSPACE].dof_count()

        globdat[gn.TIMESTEP] = self._step
        print("Running time step", self._step)

        if self._step == 0:
            globdat[gn.STATE0] = np.zeros(dc)
            globdat[gn.OLDSTATE0] = np.zeros(dc)
            globdat[gn.LAMBDA] = 0.0
            self._duOld = np.zeros(dc)

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
        for model in self.get_relevant_models("ADVANCE", self._models):
            model.ADVANCE(params, globdat)

        # Assemble K
        for model in self.get_relevant_models("GETMATRIX0", self._models):
            model.GETMATRIX0(params, globdat)

        # Assemble fext0
        fext0 = fhat * globdat[gn.LAMBDA] + self._fext0

        # Get constraints
        for model in self.get_relevant_models("GETCONSTRAINTS", self._models):
            model.GETCONSTRAINTS(params, globdat)
        cdofs, cvals = c.get_constraints()
        fdofs = [i for i in range(dc) if i not in cdofs]

        # Solve system
        duII = solveSys(K, self._fhat, c)

        if self._step < 1:
            lsign = 1
        else:  # Choose solution with smallest angle wrt previous solution
            lsign = np.sign(np.dot(self._duOld, fhat) * np.dot(duII, fhat))

        beta2F = self._beta**2 * np.dot(fhat, fhat)
        Dlam = lsign * self._dl / np.sqrt(np.dot(duII, duII) + beta2F)
        lbeta2F = Dlam * beta2F

        Du0 = Dlam * duII
        globdat[gn.STATE0] += Du0

        # Reference values to check convergence
        ref = max(Dlam * np.linalg.norm(fext0), 1)
        rel = 1.0

        # Initialize iteration loop
        while rel > self._tolerance and iteration < self._itermax:
            iteration += 1

            params[pn.MATRIX0] = np.zeros((dc, dc))
            for model in self.get_relevant_models("GETMATRIX0", self._models):
                model.GETMATRIX0(params, globdat)

            params[pn.INTFORCE] = np.zeros(dc)
            for model in self.get_relevant_models("GETINTFORCE", self._models):
                model.GETINTFORCE(params, globdat)

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
            print("Iteration %i, relative residual norm: %.4e" % (iteration, rel))

        # Alert if not convergence
        if rel > self._tolerance:
            if rel > 1:
                raise RuntimeError("Divergence in time step %i" % self._step)
            else:
                warnings.warn("No convergence in time step %i" % self._step)

        globdat[gn.LAMBDA] += Dlam

        # Check commit
        params[pn.EXTFORCE] = self._fext0
        for model in self.get_relevant_models("CHECKCOMMIT", self._models):
            model.CHECKCOMMIT(params, globdat)
        self._fext0 = params[pn.EXTFORCE]

        # Only move to next time step if commit is accepted
        if globdat[gn.ACCEPTED]:
            self._step += 1
            self._duOld = globdat[gn.STATE0] - globdat[gn.OLDSTATE0]
            globdat[gn.OLDSTATE0] = np.copy(globdat[gn.STATE0])

        while len(self._fhat) < globdat[gn.DOFSPACE].dof_count():
            self._fhat = np.append(self._fhat, 0)
            self._fext0 = np.append(self._fext0, 0)
            self._duOld = np.append(self._duOld, 0)

        if self._step >= self._nsteps:
            return "exit"
        else:
            return "ok"

    def shutdown(self, globdat):
        pass


def solveSys(K, f, c):
    Kc, fc = c.constrain(K, f)
    smat = sparse.csr_matrix(Kc)
    u = linalg.spsolve(smat, fc)
    return u

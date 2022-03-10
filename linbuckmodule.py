import numpy as np
import scipy.sparse as sparse
import scipy.sparse.linalg as linalg

from names import GlobNames as gn
from names import ParamNames as pn
from names import Actions as act

from module import *
from constrainer import Constrainer

class LinBuckModule(Module):

    def init(self, props, globdat):
        myprops = props[self._name]

    def run(self, globdat):
        dc = globdat[gn.DOFSPACE].dof_count()
        model = globdat[gn.MODEL]

        print('LinBuckModule: running unit load analysis...')

        K = np.zeros((dc, dc))
        f = np.zeros(dc)
        c = Constrainer()

        params = {}
        params[pn.MATRIX0] = K

        params[pn.CONSTRAINTS] = c

        model.take_action(act.GETMATRIX0, params, globdat)

        model.take_action(act.GETEXTFORCE, params, globdat)

        model.take_action(act.GETCONSTRAINTS, params, globdat)

        Kc, fc = c.constrain(K, f)

        smat = sparse.csr_matrix(Kc)
        u = linalg.spsolve(smat, fc)

        globdat[gn.STATE0] = u

        print('LinBuckModule: running eigenvalue problem...')


        KM = np.zeros((dc, dc))
        KG = np.zeros((dc, dc))

        params[pn.MATRIX0] = KM
        params[pn.MATRIX1] = KG

        model.take_action(act.GETMATRIXLB, params, globdat)

        cdofs,cvals = c.get_constraints()
        fdofs = [i for i in range(len(K)) if i not in cdofs]
        assert(max(max(cvals),-min(cvals))>1.e-10,'LinBuckModule does not work with nonzero Dirichlet BCs')

        lambdas,vs = eigh(KM[np._ix(cdofs,cdofs)],KG[np._ix(cdofs,cdofs)])

        for idx in np.argsort(cdofs):
            vs = np.insert(vs, cdofs[idx], cvals[idx], axis=1)

        return 'exit'

    def shutdown(self, globdat):
        pass

    def __solve(self, globdat):
        pass


def declare(factory):
    factory.declare_module('Solver', SolverModule)

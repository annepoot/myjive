import numpy as np
import scipy.sparse as sparse
import scipy.sparse.linalg as linalg
import scipy

from jive.names import GlobNames as gn
from jive.names import ParamNames as pn
from jive.names import Actions as act
from jive.app import Module
from jive.solver import Constrainer

__all__ = ["LinBuckModule"]


class LinBuckModule(Module):
    def init(self, props, globdat):
        myprops = props[self._name]

    def run(self, globdat):
        dc = globdat[gn.DOFSPACE].dof_count()
        models = globdat[gn.MODELS]

        print("LinBuckModule: running unit load analysis...")

        K = np.zeros((dc, dc))
        f = np.zeros(dc)
        globdat[gn.STATE0] = np.zeros(dc)
        c = Constrainer()

        params = {pn.MATRIX0: K, pn.CONSTRAINTS: c, pn.EXTFORCE: f}

        for model in self.get_relevant_models("GETMATRIX0", models):
            model.GETMATRIX0(params, globdat)

        for model in self.get_relevant_models("GETEXTFORCE", models):
            model.GETEXTFORCE(params, globdat)

        for model in self.get_relevant_models("GETCONSTRAINTS", models):
            model.GETCONSTRAINTS(params, globdat)

        Kc, fc = c.constrain(K, f)

        smat = sparse.csr_matrix(Kc)
        u = linalg.spsolve(smat, fc)

        globdat[gn.STATE0] = u

        print("LinBuckModule: running eigenvalue problem...")

        KM = np.zeros((dc, dc))
        KG = np.zeros((dc, dc))

        params[pn.MATRIX0] = KM
        params[pn.MATRIX1] = KG

        for model in self.get_relevant_models("GETMATRIXLB", models):
            model.GETMATRIXLB(params, globdat)

        cdofs, cvals = c.get_constraints()
        fdofs = [i for i in range(dc) if i not in cdofs]
        assert (
            max(max(cvals), -min(cvals)) < 1.0e-10
        ), "LinBuckModule does not work with nonzero Dirichlet BCs"

        ls, vs = scipy.linalg.eig(KM[np.ix_(fdofs, fdofs)], -KG[np.ix_(fdofs, fdofs)])

        for idx in np.argsort(cdofs):
            vs = np.insert(vs, cdofs[idx], cvals[idx], axis=0)

        z = list(zip(ls, vs.transpose()))
        zs = sorted(z, key=lambda f: abs(f[0]))
        lss, vss = list(zip(*zs))

        globdat[gn.LBFACTORS] = lss
        globdat[gn.HISTORY] = np.asarray(vss)
        print("LinBuckModule: critical load factor:  %8.3e" % np.real_if_close(lss[0]))

        return "exit"

    def shutdown(self, globdat):
        pass

    def __solve(self, globdat):
        pass

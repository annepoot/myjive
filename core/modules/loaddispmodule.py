import numpy as np

from jive.app import Module
from jive.names import GlobNames as gn
from jive.names import ParamNames as pn
import jive.util.proputils as pu

GROUPS = "groups"
DISP = "disp"
LOAD = "load"

__all__ = ["LoadDispModule"]


class LoadDispModule(Module):
    def init(self, props, globdat):
        if self._name in props:
            myprops = props.get(self._name)

            if GROUPS in myprops:
                self._groups = pu.parse_list(myprops[GROUPS])

        mydata = {}

        for group in self._groups:
            mydata[group] = {}
            mydata[group][LOAD] = {}
            mydata[group][DISP] = {}
            for typ in globdat[gn.DOFSPACE].get_types():
                mydata[group][LOAD][typ] = []
                mydata[group][DISP][typ] = []

        globdat[self._name] = mydata

    def run(self, globdat):
        models = globdat[gn.MODELS]
        nodes = globdat[gn.NSET]
        disp = globdat[gn.STATE0]
        dofs = globdat[gn.DOFSPACE]
        types = dofs.get_types()
        dc = dofs.dof_count()
        fint = np.zeros(dc)
        mydata = globdat[self._name]

        params = {}
        params[pn.INTFORCE] = fint

        for model in self.get_relevant_models("GETINTFORCE", models):
            model.GETINTFORCE(params, globdat)

        for group in self._groups:
            for typ in types:
                idofs = dofs.get_dofs(globdat[gn.NGROUPS][group], [typ])
                mydata[group][DISP][typ].append(np.mean(disp[idofs]))
                mydata[group][LOAD][typ].append(np.sum(fint[idofs]))

        return "ok"

    def shutdown(self, globdat):
        pass

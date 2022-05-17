import numpy as np
import matplotlib.pyplot as plt

from module import *
from names import GlobNames as gn
from names import ParamNames as pn
from names import Actions as act

import proputils as pu

XDATA = 'xData'
YDATA = 'yData'

class GraphModule(Module):

    def init(self, props, globdat):

        if self._name in props:
            myprops = props.get(self._name)

            self._xdata = pu.parse_list(myprops[XDATA])
            self._ydata = pu.parse_list(myprops[YDATA])

    def run(self, globdat):
        return 'ok'
                
    def shutdown(self, globdat):
        self._make_graph(globdat)

    def _make_graph(self, globdat):
        fig = plt.figure(2)
        for xdat, ydat in zip(self._xdata, self._ydata):
            [module, group, ld, typ] = xdat.split('.')
            x = globdat[module][group][ld][typ]
            [module, group, ld, typ] = ydat.split('.')
            y = globdat[module][group][ld][typ]
            plt.plot(x, y)
        plt.show()

def declare(factory):
    factory.declare_module('Graph', GraphModule)

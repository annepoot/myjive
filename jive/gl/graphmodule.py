import matplotlib.pyplot as plt

from jive.app.module import Module
import jive.util.proputils as pu

XDATA = "xData"
YDATA = "yData"

__all__ = ["GraphModule"]


class GraphModule(Module):
    def init(self, props, globdat):
        if self._name in props:
            myprops = props.get(self._name)

            self._xdata = pu.parse_list(myprops[XDATA])
            self._ydata = pu.parse_list(myprops[YDATA])

    def run(self, globdat):
        return "ok"

    def shutdown(self, globdat):
        self._make_graph(globdat)

    def _make_graph(self, globdat):
        fig = plt.figure(2)
        for xdat, ydat in zip(self._xdata, self._ydata):
            [module, group, ld, typ] = xdat.split(".")
            x = globdat[module][group][ld][typ]
            [module, group, ld, typ] = ydat.split(".")
            y = globdat[module][group][ld][typ]
            plt.plot(x, y)
        plt.show()

from .module import Module
from ..names import GlobNames as gn

__all__ = ["OutputModule"]


class OutputModule(Module):
    def init(self, globdat, **props):
        pass

    def run(self, globdat):
        # Temporary strat
        fname = "step" + str(globdat[gn.TIMESTEP]) + ".disp"
        u = globdat[gn.STATE0]
        print(u)

        with open(fname, "w") as out:
            for val in u:
                out.write(str(val) + "\n")

        return "ok"

    def shutdown(self, globdat):
        pass

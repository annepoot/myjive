from ..app import Module

__all__ = ["SolverModule"]


class SolverModule(Module):
    def init(self, globdat, **props):
        pass

    def run(self, globdat):
        while True:
            self.advance(globdat)

            self.solve(globdat)

            if self.commit(globdat):
                break

            print("Solution rejected; retrying\n")

            self.cancel(globdat)

        return "ok"

    def advance(self, globat):
        pass

    def solve(self, globdat):
        pass

    def cancel(self, globdat):
        pass

    def commit(self, globdat):
        return True

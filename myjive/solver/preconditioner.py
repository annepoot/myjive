from myjive.util.proputils import optional_argument

NOTIMPLEMENTEDMSG = "this function needs to be implemented in an derived class"

__all__ = ["Preconditioner", "PreconFactory"]


class PreconFactory:
    def __init__(self):
        self._precons = {}

    def declare_precon(self, typ, precon):
        self._precons[typ] = precon

    def get_precon(self, typ):
        precon = self._precons.get(typ)
        if not precon:
            raise ValueError(typ)
        return precon()

    def is_precon(self, typ):
        return typ in self._precons


class Preconditioner:

    def configure(self, globdat, **props):
        self._precision = optional_argument(self, props, "precision", default=1e-8)

    @classmethod
    def declare(cls, factory):
        name = cls.__name__
        if len(name) > 6 and name[-6:] == "Precon":
            name = name[:-6]
        factory.declare_precon(name, cls)

    def update(self, sourcematrix):
        raise NotImplementedError(NOTIMPLEMENTEDMSG)

    def dot(self, lhs):
        raise NotImplementedError(NOTIMPLEMENTEDMSG)

    def solve(self, rhs):
        raise NotImplementedError(NOTIMPLEMENTEDMSG)

    def get_matrix(self):
        raise NotImplementedError(NOTIMPLEMENTEDMSG)

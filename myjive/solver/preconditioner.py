from ..util import saveconfig as sg

from myjive.util.proputils import optional_argument

NOTIMPLEMENTEDMSG = "this function needs to be implemented in an derived class"

__all__ = ["Preconditioner", "PreconFactory"]


class PreconFactory:
    def __init__(self):
        self._precons = {}

    def declare_precon(self, typ, precon):
        self._precons[typ] = precon

    def get_precon(self, typ, name):
        precon = self._precons.get(typ)
        if not precon:
            raise ValueError(typ)
        return precon(name)

    def is_precon(self, typ):
        return typ in self._precons


class Preconditioner:
    def __init__(self, name):
        self._name = name
        self._config = {}

    @classmethod
    def declare(cls, factory):
        name = cls.__name__
        if len(name) > 6 and name[-6:] == "Precon":
            name = name[:-6]
        factory.declare_precon(name, cls)

    @sg.save_config
    def configure(self, globdat, **props):
        self._precision = optional_argument(self, props, "precision", default=1e-8)

    def update(self, sourcematrix):
        raise NotImplementedError(NOTIMPLEMENTEDMSG)

    def dot(self, lhs):
        raise NotImplementedError(NOTIMPLEMENTEDMSG)

    def solve(self, rhs):
        raise NotImplementedError(NOTIMPLEMENTEDMSG)

    def get_matrix(self):
        raise NotImplementedError(NOTIMPLEMENTEDMSG)

    def get_config(self):
        if len(self._config) == 0:
            raise NotImplementedError("Empty preconditioner get_config")
        else:
            return self._config

    def get_name(self):
        return self._name

    @staticmethod
    def save_config(configure):
        return sg.save_config(configure)

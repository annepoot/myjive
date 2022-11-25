PRECISION = 'precision'
NOTIMPLEMENTEDMSG = 'this function needs to be implemented in an derived class'

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

    def __init__(self):
        self._maxiter = 1000

    def configure(self, props):
        self._precision = props.get(PRECISION, self._precision)

    def update(self, sourcematrix):
        raise NotImplementedError(NOTIMPLEMENTEDMSG)

    def dot(self, lhs):
        raise NotImplementedError(NOTIMPLEMENTEDMSG)

    def solve(self, rhs):
        raise NotImplementedError(NOTIMPLEMENTEDMSG)

    def get_matrix(self):
        raise NotImplementedError(NOTIMPLEMENTEDMSG)

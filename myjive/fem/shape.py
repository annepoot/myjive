import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import fsolve
from warnings import warn
from numba import njit

NOTIMPLEMENTEDMSG = "this function needs to be implemented in an derived class"

__all__ = ["Shape", "ShapeFactory"]


class ShapeFactory:
    def __init__(self):
        self._creators = {}

    def declare_shape(self, typ, creator):
        self._creators[typ] = creator

    def get_shape(self, typ, ischeme):
        creator = self._creators.get(typ)
        if not creator:
            raise ValueError(typ)
        return creator(ischeme)


class Shape:
    def __init__(self, intscheme):
        # Note: these two parameters need to be implemented in the derived class
        # self._ncount = None
        # self._rank = None

        if "*" in intscheme:
            ipcount = 1
            for ip in intscheme.split("*"):
                ipcount *= int(ip.lstrip("Gauss"))
            self._int = "Gauss" + str(ipcount)
        else:
            self._int = intscheme

        if self._int.lstrip("Gauss").isnumeric():
            self._ipcount = int(self._int.lstrip("Gauss"))
        else:
            raise ValueError(self._int)

        self._ips = np.zeros((self._ipcount, self._rank))
        self._wts = np.zeros(self._ipcount)

        if self._rank == 1:
            self._ips[:, 0], self._wts[:] = leggauss(self._ipcount)
        elif self._rank == 2:
            if self._ncount == 3 or self._ncount == 6:
                if self._int == "Gauss1":
                    self._ips[0, 0] = 1.0 / 3.0
                    self._ips[0, 1] = 1.0 / 3.0
                    self._wts[0] = 0.5
                elif self._int == "Gauss3":
                    self._ips[0, 0] = 1.0 / 6.0
                    self._ips[0, 1] = 1.0 / 6.0
                    self._ips[1, 0] = 2.0 / 3.0
                    self._ips[1, 1] = 1.0 / 6.0
                    self._ips[2, 0] = 1.0 / 6.0
                    self._ips[2, 1] = 2.0 / 3.0
                    self._wts[0] = 1.0 / 6.0
                    self._wts[1] = 1.0 / 6.0
                    self._wts[2] = 1.0 / 6.0
                else:
                    raise ValueError(self._int)

            elif self._ncount == 4 or self._ncount == 9:
                sqrtipcount = np.sqrt(self._ipcount)
                if sqrtipcount.is_integer():
                    xipcount = int(sqrtipcount)
                    yipcount = int(sqrtipcount)
                else:
                    raise ValueError(self._ipcount)

                ip = 0

                for yip, ywt in zip(*leggauss(yipcount)):
                    for xip, xwt in zip(*leggauss(xipcount)):
                        self._ips[ip, 0] = xip
                        self._ips[ip, 1] = yip
                        self._wts[ip] = xwt * ywt
                        ip = ip + 1

                assert ip == self._ipcount

            else:
                raise ValueError(self._ncount)

        self._N = np.zeros((self._ipcount, self._ncount))
        self._dN = np.zeros((self._ipcount, self._rank, self._ncount))

        for ip in range(self._ipcount):
            self._N[ip] = self.eval_shape_functions(self._ips[ip])
            self._dN[ip] = self.eval_shape_gradients(self._ips[ip])

    @classmethod
    def declare(cls, factory):
        name = cls.__name__
        if len(name) > 5 and name[-5:] == "Shape":
            name = name[:-5]
        factory.declare_shape(name, cls)

    def global_rank(self):
        return self._rank

    def node_count(self):
        return self._ncount

    def ipoint_count(self):
        return self._ipcount

    def get_local_node_coords(self):
        raise NotImplementedError(NOTIMPLEMENTEDMSG)

    def get_integration_points(self):
        return self._ips

    def get_global_integration_points(self, glob_coords):
        glob_ips = np.zeros((self._ipcount, self._rank))

        for ip in range(self._ipcount):
            glob_ips[ip] = self.get_global_point(self._ips[ip], glob_coords)

        return glob_ips

    def get_integration_weights(self, glob_coords):
        wts = np.copy(self._wts)

        for ip in range(self._ipcount):
            J = self._dN[ip] @ glob_coords
            wts[ip] *= np.linalg.det(J)

        return wts

    def get_shape_functions(self):
        return self._N

    def eval_shape_functions(self, loc_point):
        raise NotImplementedError(NOTIMPLEMENTEDMSG)

    def eval_global_shape_functions(self, glob_point, glob_coords):
        loc_point = self.get_local_point(glob_point, glob_coords)
        return self.eval_shape_functions(loc_point)

    def get_global_point(self, loc_point, glob_coords):
        sfuncs = self.eval_shape_functions(loc_point)
        return sfuncs @ glob_coords

    def get_local_point(self, glob_point, glob_coords):
        # Note: since this is (in general) a non-linear problem, a non-linear solver must be called.
        # Inherited classes are encouraged to get more efficient implementations
        def f(x):
            return self.get_global_point(x, glob_coords) - glob_point

        # The initial guess is the local coordinate in the middle of the element
        x0 = np.mean(self.get_local_node_coords(), axis=0)

        # Raise an error that scipy.optimize.fsolve is necessary
        warn(
            "get_local_points needs to do a scipy.optimize.fsolve call to get a result"
        )

        # Do a non-linear solve to find the corresponding local point
        loc_point = fsolve(f, x0)

        # Make sure that the solution is actually inside the element
        if not self.contains_local_point(loc_point, tol=1e-8):
            raise ValueError(glob_point)

        return loc_point

    def contains_local_point(self, loc_point):
        raise NotImplementedError(NOTIMPLEMENTEDMSG)

    def contains_global_point(self, glob_point, glob_coords, tol=0.0):
        raise NotImplementedError(NOTIMPLEMENTEDMSG)

    def get_shape_gradients(self, glob_coords):
        return self._get_shape_gradients_jit(
            glob_coords, self._dN, self._wts, self._ipcount
        )

    @staticmethod
    @njit
    def _get_shape_gradients_jit(glob_coords, _dN, _wts, _ipcount):
        wts = np.copy(_wts)
        dN = np.copy(_dN)

        for ip in range(_ipcount):
            dNip = dN[ip]
            J = dNip @ glob_coords
            invJ, detJ = invdet(J)
            wts[ip] *= detJ
            dN[ip] = invJ @ dNip

        return dN, wts

    def eval_shape_gradients(self, loc_point):
        raise NotImplementedError(NOTIMPLEMENTEDMSG)

    def eval_global_shape_gradients(self, glob_point, glob_coords):
        loc_point = self.get_local_point(glob_point, glob_coords)
        loc_grads = self.eval_shape_gradients(loc_point)
        J = loc_grads @ glob_coords
        J_inv = np.linalg.inv(J)
        return J_inv @ loc_grads


##########################
# numba helper functions #
##########################
@njit
def det(A):
    if len(A) == 2:
        return det2x2(A)
    elif len(A) == 3:
        return det3x3(A)
    else:
        return np.linalg.det(A)


@njit
def det2x2(A):
    return det2x2i(A[0, 0], A[0, 1], A[1, 0], A[1, 1])


@njit
def det2x2i(a, b, c, d):
    return a * d - b * c


@njit
def det3x3(A):
    return (
        A[0, 0] * A[1, 1] * A[2, 2]
        + A[0, 1] * A[1, 2] * A[2, 0]
        + A[0, 2] * A[1, 0] * A[2, 1]
        - A[0, 0] * A[1, 2] * A[2, 1]
        - A[0, 1] * A[1, 0] * A[2, 2]
        - A[0, 2] * A[1, 1] * A[2, 0]
    )


@njit
def invdet(A):
    if len(A) == 1:
        return invdet1x1(A)
    elif len(A) == 2:
        return invdet2x2(A)
    elif len(A) == 3:
        return invdet3x3(A)
    else:
        raise ValueError("array is too big for invdet function")


@njit
def invdet1x1(A):
    a = A[0, 0]
    return np.array([[1.0 / a]]), a


@njit
def invdet2x2(A):
    a, b, c, d = A[0, 0], A[0, 1], A[1, 0], A[1, 1]
    det = det2x2i(a, b, c, d)
    adj = np.array([[d, -b], [-c, a]])
    return adj / det, det


@njit
def cof3x3(A, r, c):
    return det2x2i(
        A[(r - 1) % 3, (c - 1) % 3],
        A[(r - 1) % 3, (c + 1) % 3],
        A[(r + 1) % 3, (c - 1) % 3],
        A[(r + 1) % 3, (c + 1) % 3],
    )


@njit
def invdet3x3(A):
    det = det3x3(A)
    adj = np.array(
        [
            [cof3x3(A, 0, 0), cof3x3(A, 1, 0), cof3x3(A, 2, 0)],
            [cof3x3(A, 0, 1), cof3x3(A, 1, 1), cof3x3(A, 2, 1)],
            [cof3x3(A, 0, 2), cof3x3(A, 1, 2), cof3x3(A, 2, 2)],
        ]
    )
    return adj / det, det

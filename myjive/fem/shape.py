import numpy as np
from scipy.optimize import fsolve
from warnings import warn
from .jit.shape import get_shape_gradients_jit

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

        self._int = intscheme

        if self._int == "Gauss1":
            self._ipcount = 1
        elif self._int == "Gauss2":
            self._ipcount = 2
        elif self._int == "Gauss3":
            self._ipcount = 3
        elif self._int == "Gauss4":
            self._ipcount = 4
        elif self._int == "Gauss9":
            self._ipcount = 9
        else:
            raise ValueError(self._int)

        self._ips = np.zeros((self._rank, self._ipcount))
        self._wts = np.zeros(self._ipcount)

        if self._rank == 1:
            if self._int == "Gauss1":
                self._ips[0, 0] = 0
                self._wts[0] = 2

            elif self._int == "Gauss2":
                self._ips[0, 0] = -1 / np.sqrt(3)
                self._ips[0, 1] = 1 / np.sqrt(3)
                self._wts[0] = 1
                self._wts[1] = 1

            elif self._int == "Gauss3":
                self._ips[0, 0] = -np.sqrt(3.0 / 5.0)
                self._ips[0, 1] = 0
                self._ips[0, 2] = np.sqrt(3.0 / 5.0)
                self._wts[0] = 5.0 / 9.0
                self._wts[1] = 8.0 / 9.0
                self._wts[2] = 5.0 / 9.0

            else:
                raise ValueError(self._int)

        elif self._rank == 2:
            if self._ncount == 3 or self._ncount == 6:
                if self._int == "Gauss1":
                    self._ips[0, 0] = 1.0 / 3.0
                    self._ips[1, 0] = 1.0 / 3.0
                    self._wts[0] = 0.5
                elif self._int == "Gauss3":
                    self._ips[0, 0] = 1.0 / 6.0
                    self._ips[1, 0] = 1.0 / 6.0
                    self._ips[0, 1] = 2.0 / 3.0
                    self._ips[1, 1] = 1.0 / 6.0
                    self._ips[0, 2] = 1.0 / 6.0
                    self._ips[1, 2] = 2.0 / 3.0
                    self._wts[0] = 1.0 / 6.0
                    self._wts[1] = 1.0 / 6.0
                    self._wts[2] = 1.0 / 6.0

                else:
                    raise ValueError(self._int)

            elif self._ncount == 4 or self._ncount == 9:
                if self._int == "Gauss1":
                    self._ips[0, 0] = 0.0
                    self._ips[1, 0] = 0.0
                    self._wts[0] = 4.0
                elif self._int == "Gauss4":
                    invsqrt3 = 1 / np.sqrt(3)
                    self._ips[0, 0] = -invsqrt3
                    self._ips[1, 0] = -invsqrt3
                    self._ips[0, 1] = invsqrt3
                    self._ips[1, 1] = -invsqrt3
                    self._ips[0, 2] = invsqrt3
                    self._ips[1, 2] = invsqrt3
                    self._ips[0, 3] = -invsqrt3
                    self._ips[1, 3] = invsqrt3

                    self._wts[0] = 1.0
                    self._wts[1] = 1.0
                    self._wts[2] = 1.0
                    self._wts[3] = 1.0
                elif self._int == "Gauss9":
                    invsqrt35 = 1 / np.sqrt(3.0 / 5.0)
                    self._ips[0, 0] = -invsqrt35
                    self._ips[1, 0] = -invsqrt35
                    self._ips[0, 1] = 0
                    self._ips[1, 1] = -invsqrt35
                    self._ips[0, 2] = invsqrt35
                    self._ips[1, 2] = -invsqrt35
                    self._ips[0, 3] = -invsqrt35
                    self._ips[1, 3] = 0
                    self._ips[0, 4] = 0
                    self._ips[1, 4] = 0
                    self._ips[0, 5] = invsqrt35
                    self._ips[1, 5] = 0
                    self._ips[0, 6] = -invsqrt35
                    self._ips[1, 6] = invsqrt35
                    self._ips[0, 7] = 0
                    self._ips[1, 7] = invsqrt35
                    self._ips[0, 8] = invsqrt35
                    self._ips[1, 8] = invsqrt35

                    self._wts[0] = 25.0 / 81.0
                    self._wts[1] = 40.0 / 81.0
                    self._wts[2] = 25.0 / 81.0
                    self._wts[3] = 40.0 / 81.0
                    self._wts[4] = 64.0 / 81.0
                    self._wts[5] = 40.0 / 81.0
                    self._wts[6] = 25.0 / 81.0
                    self._wts[7] = 40.0 / 81.0
                    self._wts[8] = 25.0 / 81.0

                else:
                    raise ValueError(self._int)

            else:
                raise ValueError(self._ncount)

        self._N = np.zeros((self._ncount, self._ipcount))
        self._dN = np.zeros((self._ncount, self._rank, self._ipcount))

        for ip in range(self._ipcount):
            self._N[:, ip] = self.eval_shape_functions(self._ips[:, ip])
            self._dN[:, :, ip] = self.eval_shape_gradients(self._ips[:, ip])

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
        glob_ips = np.zeros((self._rank, self._ipcount))

        for ip in range(self._ipcount):
            glob_ips[:, ip] = self.get_global_point(self._ips[:, ip], glob_coords)

        return glob_ips

    def get_integration_weights(self, glob_coords):
        wts = np.copy(self._wts)

        for ip in range(self._ipcount):
            J = np.matmul(glob_coords, self._dN[:, :, ip])
            wts[ip] *= np.linalg.det(J)

        return wts

    def get_shape_functions(self):
        return self._N

    def eval_shape_functions(self, loc_point):
        raise NotImplementedError(NOTIMPLEMENTEDMSG)

    def get_global_point(self, loc_point, glob_coords):
        sfuncs = self.eval_shape_functions(loc_point)
        return np.matmul(glob_coords, sfuncs)

    def get_local_point(self, glob_point, glob_coords):
        # Note: since this is (in general) a non-linear problem, a non-linear solver must be called.
        # Inherited classes are encouraged to get more efficient implementations
        def f(x):
            return self.get_global_point(x, glob_coords) - glob_point

        # The initial guess is the local coordinate in the middle of the element
        x0 = np.mean(self.get_local_node_coords(), axis=1)

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

    def get_shape_gradients(self, glob_coords):
        return get_shape_gradients_jit(glob_coords, self._dN, self._wts, self._ipcount)

    def eval_shape_gradients(self, loc_point):
        raise NotImplementedError(NOTIMPLEMENTEDMSG)

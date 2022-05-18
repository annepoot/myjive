
import numpy as np

NOTIMPLEMENTEDMSG = 'this function needs to be implemented in an derived class'

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

        if self._int == 'Gauss1':
            self._ipcount = 1
        elif self._int == 'Gauss2':
            self._ipcount = 2

        self._ips = np.zeros((self._rank, self._ipcount))
        self._wts = np.zeros(self._ipcount)

        if self._rank == 1:

            if self._int == 'Gauss1':
                self._ips[0, 0] = 0
                self._wts[0] = 2

            elif self._int == 'Gauss2':
                self._ips[0, 0] = -1 / np.sqrt(3)
                self._ips[0, 1] = 1 / np.sqrt(3)
                self._wts[0] = 1
                self._wts[1] = 1

        elif self._rank == 2:
            if self._int == 'Gauss1':
                self._ips[0, 0] = 1.0 / 3.0
                self._ips[1, 0] = 1.0 / 3.0
                self._wts[0] = 0.5

        self._N = np.zeros((self._ncount, self._ipcount))
        self._dN = np.zeros((self._ncount, self._rank, self._ipcount))

        for ip in range(self._ipcount):
            self._N[:,ip] = self.eval_shape_functions(self._ips[:,ip])
            self._dN[:,:,ip] = self.eval_shape_gradients(self._ips[:,ip])


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
            glob_ips[:,ip] = self.get_global_point(self._ips[:,ip], glob_coords)

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

    def get_local_point(self):
        raise NotImplementedError(NOTIMPLEMENTEDMSG)

    def contains_local_point(self):
        raise NotImplementedError(NOTIMPLEMENTEDMSG)

    def get_shape_gradients(self, glob_coords):
        wts = np.copy(self._wts)
        dN = np.copy(self._dN)

        for ip in range(self._ipcount):
            J = np.matmul(glob_coords, dN[:, :, ip])
            wts[ip] *= np.linalg.det(J)
            dN[:, :, ip] = np.matmul(dN[:, :, ip], np.linalg.inv(J))

        return dN, wts

    def eval_shape_gradients(self, loc_point):
        raise NotImplementedError(NOTIMPLEMENTEDMSG)

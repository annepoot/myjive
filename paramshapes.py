import numpy as np

from shape import *


class Tri3Shape(Shape):
    def __init__(self, intscheme):
        print('Creating Tri3Shape...')

        self._ncount = 3
        self._rank = 2
        self._int = intscheme

        if self._int == 'Gauss1':
            self._ipcount = 1
            self._ips = np.zeros((2, 1))
            self._wts = np.zeros(1)

            self._ips[0, 0] = 1.0 / 3.0
            self._ips[1, 0] = 1.0 / 3.0
            self._wts[0] = 0.5
        else:
            raise SyntaxError('Unsupported integration scheme for Triangle3 shape')

        # TODO: remove middle index from N?

        self._N = np.zeros((ncount, 1, self._ipcount))
        self._dN = np.zeros((ncount, rank, self._ipcount))

        self._N[0, :, :] = 1.0 - self._ips[0, :] - self._ips[1, :]
        self._N[1, :, :] = self._ips[0, :]
        self._N[2, :, :] = self._ips[1, :]

        self._dN[0, 0, :] = -1.0
        self._dN[0, 1, :] = -1.0
        self._dN[1, 0, :] = 1.0
        self._dN[1, 1, :] = 0.0
        self._dN[2, 0, :] = 0.0
        self._dN[2, 1, :] = 1.0

    def get_shape_gradients(self, coords):
        wts = self._wts
        dN = self._dN

        for ip in range(self._ipcount):
            J = np.matmul(coords, dN[:, :, ip])
            wts[ip] *= np.linalg.det(J)
            dN[:, :, ip] = np.matmul(dN[:, :, ip], np.linalg.inv(J))

        return dN, wts

    def get_shape_functions(self):
        return self._N


class Line2Shape(Shape):
    def __init__(self, intscheme):
        print('Creating Line2Shape...')

        self._ncount = 2
        self._rank = 1
        self._int = intscheme

        if self._int == 'Gauss1':
            self._ipcount = 1
            self._ips = np.zeros((1, 1))
            self._wts = np.zeros(1)

            self._ips[0, 0] = 0
            self._wts[0] = 2
        else:
            raise SyntaxError('Unsupported integration scheme for Line2 shape')

        self._N = np.zeros((2, 1, self._ipcount))
        self._dN = np.zeros((2, 1, self._ipcount))

        self._N[0, :, :] = 0.5 - 0.5 * self._ips[0, :]
        self._N[1, :, :] = 0.5 + 0.5 * self._ips[0, :]

        self._dN[0, 0, :] = -0.5
        self._dN[1, 0, :] = 0.5

    def get_shape_gradients(self, coords):
        wts = self._wts
        dN = self._dN

        for ip in range(self._ipcount):
            J = np.matmul(coords, dN[:, :, ip]).item()
            wts[ip] *= J
            dN[:, :, ip] = dN[:, :, ip] / J

        return dN, wts

    def get_shape_functions(self):
        return self._N


def declare(factory):
    factory.declare_shape('Triangle3', Tri3Shape)
    factory.declare_shape('Line2', Line2Shape)

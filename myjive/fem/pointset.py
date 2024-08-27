import numpy as np

from .itemset import ItemSet, XItemSet

__all__ = ["PointSet", "XPointSet", "to_xpointset"]


class PointSet(ItemSet):
    def __init__(self, points=None):
        super().__init__(points)
        self._rank = 0

    def rank(self):
        return self._rank

    def get_point_coords(self, ipoint):
        return self._data[ipoint]

    def get_coords(self):
        coords = []
        for ipoint in range(self.size()):
            coords.append(self.get_point_coords(ipoint))
        return np.array(coords)

    def get_some_coords(self, ipoints):
        coords = []
        for ipoint in ipoints:
            coords.append(self.get_point_coords(ipoint))
        return np.array(coords)


class XPointSet(PointSet, XItemSet):
    def add_point(self, coords):
        if self.size() == 0:
            self._rank = len(coords)
        else:
            assert self._rank == len(coords)
        self.add_item(coords)

    def erase_point(self, ipoint):
        self.erase_item(ipoint)

    def set_point_coords(self, ipoint, coords):
        self._data[ipoint] = coords

    def set_points(self, coords):
        if coords.shape[0] != self.size():
            raise ValueError(
                "first dimension of coords does not match the number of points"
            )
        for ipoint in range(self.size()):
            self.set_point_coords(ipoint, coords[ipoint])

    def set_some_coords(self, ipoints, coords):
        if coords.shape[0] != self.size():
            raise ValueError(
                "first dimension of coords does not match the size of ipoints"
            )
        for i, ipoint in enumerate(ipoints):
            self.set_point_coords(ipoint, coords[i])

    def to_pointset(self):
        self.__class__ = PointSet
        return self


def to_xpointset(points):
    points.__class__ = XPointSet
    return points

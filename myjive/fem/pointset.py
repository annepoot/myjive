import numpy as np

from .itemset import ItemSet, XItemSet, ItemMap

__all__ = ["PointSet", "XPointSet", "to_xpointset"]


class PointSet(ItemSet):
    def __init__(self, points=None):
        if points is None:
            self._data = np.zeros((0, 0))
            self._map = ItemMap()
            self._size = 0
            self._rank = 0
        else:
            self._data = points._data
            self._map = points._map
            self._size = points._size
            self._rank = points._rank

    def __len__(self):
        return self._size

    def __iter__(self):
        return iter(self._data[: self._size])

    def __next__(self):
        return next(self._data[: self._size])

    def __getitem__(self, ipoint):
        return self.get_point_coords(ipoint)

    def size(self):
        return self._size

    def rank(self):
        return self._rank

    def get_point_coords(self, ipoint):
        return self._data[ipoint]

    def get_coords(self):
        return self._data[: self._size]

    def get_some_coords(self, ipoints):
        return self._data[ipoints]


class XPointSet(PointSet, XItemSet):
    def add_point(self, coords, point_id=None):
        if self._size + 1 > len(self._data):
            if self._size == 0:
                self._rank = len(coords)
                self._data = np.zeros((1, self._rank))
            else:
                # self._data = self._data.resize((2 * self._size, self._rank))
                self._data = np.pad(self._data, ((0, self._size), (0, 0)))

        self._data[self._size] = coords
        self._size += 1

        self._map.add_item(point_id)

        return self._size - 1

    def add_points(self, coords, point_ids=None):
        assert len(coords.shape) == 2
        add_node_count = coords.shape[0]

        if self._size + add_node_count > len(self._data):
            if self._size == 0:
                self._rank = coords.shape[1]
                self._data = np.zeros((add_node_count, self._rank))
            else:
                padding = max(self._size, add_node_count)
                self._data = np.pad(self._data, ((0, padding), (0, 0)))

        self._data[self._size : self._size + add_node_count] = coords
        self._size += add_node_count

        self._map.add_items(add_node_count, point_ids)

        return np.arange(self._size - add_node_count, self._size)

    def erase_point(self, ipoint):
        self._data = np.delete(self._data, ipoint, axis=0)
        self._map.erase_item(ipoint)
        self._size -= 1

    def erase_points(self, ipoints):
        self._data = np.delete(self._data, ipoints, axis=0)
        self._map.erase_item(ipoints)
        self._size -= len(ipoints)

    def set_point_coords(self, ipoint, coords):
        self._data[ipoint] = coords

    def set_points(self, coords):
        if coords.shape[0] != self.size():
            raise ValueError(
                "first dimension of coords does not match the number of points"
            )
        self._data = coords
        self._rank = self._data.shape[1]

    def set_some_coords(self, ipoints, coords):
        raise NotImplementedError("has not been tested yet")
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

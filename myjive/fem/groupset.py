import numpy as np

from .itemset import ItemSet, XItemSet, ItemMap

__all__ = ["GroupSet", "XGroupSet", "to_xgroupset"]


class GroupSet(ItemSet):
    def __init__(self, groups=None):
        if groups is None:
            self._data = np.zeros((0, 0), dtype=int)
            self._groupsizes = np.zeros(0, dtype=int)
            self._map = ItemMap()
            self._size = 0
            self._maxgroupsize = 0
        else:
            self._data = groups._data
            self._groupsizes = groups._groupsizes
            self._map = groups._map
            self._size = groups._size
            self._maxgroupsize = groups._maxgroupsize

    def __len__(self):
        return self._size

    def __iter__(self):
        for igroup in range(self._size):
            groupsize = self._groupsizes[igroup]
            yield self._data[igroup][:groupsize]

    def __getitem__(self, igroup):
        return self.get_group_members(igroup)

    def size(self):
        return self._size

    def get_group_size(self, igroup):
        return self._groupsizes[igroup]

    def get_group_members(self, igroup):
        return self._data[igroup][: self.get_group_size(igroup)]

    def get_some_group_members(self, index, igroup):
        return self.get_group_members(igroup)[index]

    def get_members_of(self, igroups):
        members_of = np.zeros(self._maxgroupsize * len(igroups), dtype=int)
        offset = 0
        for igroup in igroups:
            members = self.get_group_members(igroup)
            groupsize = len(members)
            members_of[offset : offset + groupsize] = members
            offset += groupsize
        return members_of[:offset]

    def get_unique_members_of(self, igroups):
        return np.unique(self.get_members_of(igroups))

    def max_group_size(self):
        return self._maxgroupsize

    def max_group_size_of(self, igroups):
        return max(self._groupsizes[igroups])


class XGroupSet(GroupSet, XItemSet):
    def add_group(self, members, group_id=None):
        if self._size == len(self._data):
            if self._size == 0:
                self._data = np.zeros((1, 0), dtype=int)
                self._groupsizes = np.zeros(1, dtype=int)
            else:
                self._data = np.pad(self._data, ((0, self._size), (0, 0)))
                self._groupsizes = np.pad(self._groupsizes, (0, self._size))

        groupsize = len(members)
        if groupsize > self._maxgroupsize:
            self._data = np.pad(
                self._data, ((0, 0), (0, groupsize - self._maxgroupsize))
            )
            self._maxgroupsize = groupsize

        self._data[self._size, :groupsize] = members
        self._groupsizes[self._size] = groupsize
        self._size += 1

        self._map.add_item(group_id)

        return self._size - 1

    def erase_group(self, igroup):
        self._data = np.delete(self._data, igroup, axis=0)
        self._groupsizes = np.delete(self._groupsizes, igroup, axis=0)
        self._map.erase_item(igroup)
        self._size -= 1

    def set_group_members(self, igroup, members):
        raise NotImplementedError("has not been tested yet")
        groupsize = len(members)
        self._data[igroup, :groupsize] = members
        self._groupsizes[igroup] = groupsize

    def to_groupset(self):
        self.__class__ = GroupSet
        return self


def to_xgroupset(points):
    points.__class__ = XGroupSet
    return points

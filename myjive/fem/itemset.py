import numpy as np

__all__ = ["ItemSet", "XItemSet"]


class ItemSet:
    def __init__(self, items=None):
        if items is None:
            self._data = []
            self._map = ItemMap()
        else:
            self._data = items._data
            self._map = items._map

    def __len__(self):
        return self.size()

    def __iter__(self):
        return iter(self._data)

    def __next__(self):
        return next(self._data)

    def __getitem__(self, iitem):
        return self._data[iitem]

    def size(self):
        return len(self._data)

    def get_item_map(self):
        return self._map

    def find_item(self, item_id):
        return self._map.find_item(item_id)

    def find_items(self, item_ids):
        return self._map.find_items(item_ids)

    def get_item_id(self, iitem):
        return self._map.get_item_id(iitem)

    def get_item_ids(self, iitems):
        return self._map.get_item_ids(iitems)


class XItemSet(ItemSet):
    def clear(self):
        self._data = []
        self._map = ItemMap()

    def add_item(self, item, item_id=None):
        self._data.append(item)
        self._map.add_item(item_id)

    def erase_item(self, iitem):
        self._data.pop(iitem)
        self._map.erase_item(iitem)


class ItemMap:
    def __init__(self):
        self._map = {}
        self._invmap = []
        self._maxkey = None

    def get_item_maps(self):
        return self._map, self._invmap

    def find_item(self, item_id):
        return self._map.get(item_id, -1)

    def find_items(self, item_ids):
        iitems = np.empty_like(item_ids, dtype=int)
        for i, item_id in enumerate(item_ids):
            iitems[i] = self.find_item(item_id)
        return iitems

    def get_item_id(self, iitem):
        return self._invmap[iitem]

    def get_item_ids(self, iitems):
        return list(map(self._invmap.__getitem__, list(iitems)))

    def clear(self):
        self._map = {}
        self._invmap = []

    def add_item(self, item_id=None):
        size = len(self._map)

        if item_id is None:
            if size == 0:
                item_id = 1
            else:
                if self._maxkey is None:
                    self._maxkey = max(
                        [i for i in self._map.keys() if isinstance(i, int)]
                    )
                item_id = self._maxkey + 1
            self._maxkey = item_id
        else:
            self._maxkey = None

        if item_id in self._map.keys():
            raise ValueError("item ID already exists in itemset")

        self._map[item_id] = size
        self._invmap.append(item_id)

    def add_items(self, n_item, item_ids=None):
        size = len(self._map)

        if item_ids is None:
            if size == 0:
                item_ids = np.arange(1, n_item + 1)
            else:
                if self._maxkey is None:
                    self._maxkey = max(
                        [i for i in self._map.keys() if isinstance(i, int)]
                    )
                item_ids = np.arange(self._maxkey + 1, self._maxkey + n_item + 1)
            self._maxkey = item_ids[-1]

        newdict = dict(zip(item_ids, np.arange(size, size + n_item)))

        if len(newdict.keys() & self._map.keys()) > 0:
            raise ValueError("item ID already exists in itemset")

        self._map.update(newdict)
        self._invmap.extend(item_ids)

    def erase_item(self, iitem):
        for item_id, idx in self._map.items():
            if idx > iitem:
                self._map[item_id] = idx - 1
            elif idx == iitem:
                pop_id = item_id
        self._map.pop(pop_id)
        self._invmap.pop(iitem)

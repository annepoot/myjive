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
    def __init__(self, imap=None):
        if imap is None:
            self._map = {}
        else:
            self._map = imap

    def get_item_map(self):
        return self._map

    def find_item(self, item_id):
        return self._map.get(item_id, -1)

    def find_items(self, item_ids):
        iitems = np.empty_like(item_ids, dtype=int)
        for i, item_id in enumerate(item_ids):
            iitems[i] = self.find_item(item_id)
        return np.array(iitems, dtype=int)

    def get_item_id(self, iitem):
        for item_id, idx in self._map.items():
            if idx == iitem:
                return item_id
        else:
            return -1

    def get_item_ids(self, iitems):
        item_ids = np.empty_like(iitems, dtype=int)
        for i, iitem in enumerate(iitems):
            item_ids[i] = self.get_item_id(iitem)
        return np.array(item_ids, dtype=int)

    def clear(self):
        self._map = {}

    def add_item(self, item_id=None):
        size = len(self._map)
        if item_id is None:
            if size == 0:
                item_id = 1
            else:
                maxkey = max(self._map, key=self._map.get)
                item_id = type(maxkey)(int(maxkey) + 1)
        if item_id in self._map.keys():
            raise ValueError("item ID already exists in itemset")
        self._map[item_id] = size

    def erase_item(self, iitem):
        for item_id, idx in self._map.items():
            if idx > iitem:
                self._map[item_id] = idx - 1
            elif idx == iitem:
                pop_id = item_id
        self._map.pop(pop_id)

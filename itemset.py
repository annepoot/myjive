import numpy as np

class ItemSet():

    def __init__(self):
        self._map = {}
        self._data = []

    def size(self):
        return len(self._data)

    def get_item_map(self):
        return self._map


    def find_item(self, item_id):
        return self._map.get(item_id, -1)

    def find_items(self, item_ids):
        iitems = np.empty_like(item_ids, dtype=int)
        for i, item_id in enumerate(item_ids):
            iitems[i] = self.find_item(item_id)
        return iitems

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
        return item_ids


class XItemSet(ItemSet):

    def clear(self):
        self._map = {}
        self._data= []

    def add_item(self, item, item_id=None):
        size = self.size()
        if item_id is None:
            item_id = size+1
        if item_id in self._map.keys():
            raise ValueError('item ID already exists in itemset')
        self._map[item_id] = size+1
        self._data.append(item)

    def erase_item(self, iitem):
        for item_id, idx in self._map.items():
            if idx > iitem:
                self._map[item_id] = idx - 1
            elif idx == iitem:
                self._map.pop(item_id)
        self._data.pop(iitem)

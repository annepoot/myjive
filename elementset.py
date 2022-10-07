import numpy as np
from itemset import ItemSet, XItemSet
from element import Element

class ElementSet(ItemSet):

    def __init__(self, data=None):
        super().__init__(data)
        self._maxNodeCount = 0

    def find_element(self, elem_id):
        return self.find_item(elem_id)

    def get_elem_id(self, ielem):
        return self.get_item_id(ielem)

    def max_elem_node_count(self):
        return self._maxNodeCount

    def max_elem_node_count_of(self, ielems):
        maxNodeCount = 0
        for ielem in ielems:
            maxNodeCount = max(maxNodeCount, self.get_elem_node_count(ielem))
        return maxNodeCount

    def get_elem_node_count(self, ielem):
        return self._data[ielem].get_node_count()

    def get_elem_nodes(self, ielem):
        return self._data[ielem].get_nodes()

    def get_nodes(self):
        nodes = []
        for ielem in range(self.size()):
            nodes.append(self.get_elem_nodes(ielem))
        return nodes

    def get_some_elem_nodes(self, index, inode):
        return self.get_elem_nodes(inode)[index]


class XElementSet(ElementSet, XItemSet):

    def add_element(self, inodes, elem_id=None):
        elem = Element(inodes)
        nodeCount = elem.get_node_count()

        if self._maxNodeCount < nodeCount:
            self._maxNodeCount = nodeCount

        self.add_item(elem, elem_id)

    def erase_element(self, ielem):
        nodeCount = self.get_elem_node_count(ielem)

        self.erase_item(ielem)

        if nodeCount == self.max_elem_node_count():
            for elem in self._data:
                if elem.get_node_count() == nodeCount:
                    break
            else:
                self._maxNodeCount = self.max_elem_node_count_of(range(self.size()))

    def set_elem_nodes(self, ielem, nodes):
        self._data[ielem].set_nodes(nodes)

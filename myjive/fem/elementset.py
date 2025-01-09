from .groupset import GroupSet, XGroupSet

__all__ = ["ElementSet", "XElementSet", "to_xelementset"]


class ElementSet(GroupSet):
    def __init__(self, nodes, elems=None):
        super().__init__(elems)
        self._nodes = nodes

    def get_nodes(self):
        return self._nodes

    def find_element(self, elem_id):
        return self.find_item(elem_id)

    def get_elem_id(self, ielem):
        return self.get_item_id(ielem)

    def max_elem_node_count(self):
        return self.max_group_size()

    def max_elem_node_count_of(self, ielems):
        return self.max_group_size_of(ielems)

    def get_elem_node_count(self, ielem):
        return self.get_group_size(ielem)

    def get_elem_nodes(self, ielem):
        return self.get_group_members(ielem)

    def get_some_elem_nodes(self, index, ielem):
        return self.get_some_group_members(index, ielem)

    def get_nodes_of(self, ielems):
        return self.get_members_of(ielems)

    def get_unique_nodes_of(self, ielems):
        return self.get_unique_members_of(ielems)


class XElementSet(ElementSet, XGroupSet):
    def add_element(self, inodes, elem_id=None):
        return self.add_group(inodes, elem_id)

    def add_elements(self, inodes, sizes, elem_ids=None):
        return self.add_groups(inodes, sizes, elem_ids)

    def erase_element(self, ielem):
        self.erase_group(ielem)

    def erase_elements(self, ielems):
        self.erase_groups(ielems)

    def set_elem_nodes(self, ielem, nodes):
        self.set_group_members(ielem, nodes)

    def to_elementset(self):
        self.__class__ = ElementSet
        return self


def to_xelementset(elems):
    elems.__class__ = XElementSet
    return elems

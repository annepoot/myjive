from .itemgroup import ItemGroup

__all__ = ["ElementGroup"]


class ElementGroup(ItemGroup):
    def __init__(self, elements, data=None):
        from .elementset import ElementSet

        super().__init__(elements, data)

        assert isinstance(self._items, ElementSet)

    def get_elements(self):
        return self._items

    def get_node_indices(self):
        return self._items.get_unique_nodes_of(self.get_indices())

    def add_element(self, ielem):
        self.add_item(ielem)

    def add_elements(self, ielems):
        self.add_items(ielems)

    def erase_element(self, ielem):
        self.erase_item(ielem)

    def erase_elements(self, ielems):
        self.erase_items(ielems)

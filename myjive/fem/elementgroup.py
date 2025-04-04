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

from .itemgroup import ItemGroup

__all__ = ["NodeGroup"]


class NodeGroup(ItemGroup):
    def __init__(self, nodes, data=None):
        from .nodeset import NodeSet

        super().__init__(nodes, data)

        assert isinstance(self._items, NodeSet)

    def get_nodes(self):
        return self._items

    def get_coords(self):
        return self._items.get_some_coords(self.get_indices())

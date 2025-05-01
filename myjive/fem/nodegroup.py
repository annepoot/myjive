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

    def add_node(self, inode):
        self.add_item(inode)

    def add_nodes(self, inodes):
        self.add_items(inodes)

    def erase_node(self, inode):
        self.erase_item(inode)

    def erase_nodes(self, inodes):
        self.erase_items(inodes)

from .pointset import PointSet, XPointSet

__all__ = ["NodeSet", "XNodeSet", "to_xnodeset"]


class NodeSet(PointSet):
    def find_node(self, node_id):
        return self.find_item(node_id)

    def find_nodes(self, node_ids):
        return self.find_items(node_ids)

    def get_node_id(self, inode):
        return self.get_item_id(inode)

    def get_node_ids(self, inodes):
        return self.get_item_ids(inodes)

    def get_node_coords(self, inode):
        return self.get_point_coords(inode)


class XNodeSet(NodeSet, XPointSet):
    def add_node(self, coords, node_id=None):
        return self.add_point(coords, node_id)

    def add_nodes(self, coords, node_ids=None):
        return self.add_points(coords, node_ids)

    def erase_node(self, inode):
        self.erase_point(inode)

    def erase_nodes(self, inodes):
        self.erase_points(inodes)

    def set_node_coords(self, inode, coords):
        self.set_point_coords(inode, coords)

    def set_coords(self, coords):
        self.set_points(coords)

    def set_some_coords(self, inodes, coords):
        self.set_some_points(inodes, coords)

    def to_nodeset(self):
        self.__class__ = NodeSet
        return self


def to_xnodeset(nodes):
    nodes.__class__ = XNodeSet
    return nodes

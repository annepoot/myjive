class ShapeFactory:
    def __init__ (self):
        self._creators = {}

    def declare_shape (self, typ, creator):
        self._creators[typ] = creator

    def get_shape (self, typ, ischeme):
        creator = self._creators.get (typ)
        if not creator:
            raise ValueError (typ)
        return creator(ischeme)

class Shape:
    def __init__ (self, intscheme):
        pass

    def node_count (self):
        return self._ncount

    def ipoint_count (self):
        return self._ipcount

    def global_rank (self):
        return self._rank

    def get_shape_gradients (self, coords):
        pass


import numpy as np

from shape import Shape

class Tri3Shape(Shape):
    def __init__(self, intscheme):
        print('Creating Tri3Shape...')

        # Set the nodecount and rank of the elements
        self._ncount = 3
        self._rank = 2

        # Refer to the Shape class to handle the rest of the initialization
        super().__init__(intscheme)

    def get_local_node_coords(self):
        # Return the standard triangle with nodes at (0,0), (1,0) and (0,1)
        loc_coords = np.zeros((self._rank, self._ncount))

        loc_coords[0, 0] = 0.0
        loc_coords[0, 1] = 1.0
        loc_coords[0, 2] = 0.0
        loc_coords[1, 0] = 0.0
        loc_coords[1, 1] = 0.0
        loc_coords[1, 2] = 1.0

        return loc_coords

    def eval_shape_functions(self, loc_point):
        # Evalulate the shape functions in the local coordinate system
        sfuncs = np.zeros(self._ncount)

        sfuncs[0] = 1.0 - loc_point[0] - loc_point[1]
        sfuncs[1] = loc_point[0]
        sfuncs[2] = loc_point[1]

        return sfuncs

    def eval_shape_gradients(self, loc_point):
        # Evaluate the shape gradients in the local coordinate system
        # Note that no weights are applied!
        sgrads = np.zeros((self._ncount, self._rank))

        sgrads[0, 0] = -1.0
        sgrads[0, 1] = -1.0
        sgrads[1, 0] = 1.0
        sgrads[1, 1] = 0.0
        sgrads[2, 0] = 0.0
        sgrads[2, 1] = 1.0

        return sgrads


class Tri6Shape(Shape):
    def __init__(self, intscheme):
        print('Creating Tri6Shape...')

        # Set the nodecount and rank of the elements
        self._ncount = 6
        self._rank = 2

        # Refer to the Shape class to handle the rest of the initialization
        super().__init__(intscheme)

    def get_local_node_coords(self):
        # Return the standard triangle with nodes at (0,0), (1,0) and (0,1)
        loc_coords = np.zeros((self._rank, self._ncount))

        loc_coords[0, 0] = 0.0
        loc_coords[0, 1] = 0.5
        loc_coords[0, 2] = 1.0
        loc_coords[0, 3] = 0.5
        loc_coords[0, 4] = 0.0
        loc_coords[0, 5] = 0.0
        loc_coords[1, 0] = 0.0
        loc_coords[1, 1] = 0.0
        loc_coords[1, 2] = 0.0
        loc_coords[1, 3] = 0.5
        loc_coords[1, 4] = 1.0
        loc_coords[1, 5] = 0.5

        return loc_coords

    def eval_shape_functions(self, loc_point):
        # Evalulate the shape functions in the local coordinate system
        sfuncs = np.zeros(self._ncount)

        sfuncs[0] = 2 * (0.5 - loc_point[0] - loc_point[1]) * (1 - loc_point[0] - loc_point[1])
        sfuncs[1] = 4 * loc_point[0] * (1 - loc_point[0] - loc_point[1])
        sfuncs[2] = -2 * loc_point[0] * (0.5 - loc_point[0])
        sfuncs[3] = 4 * loc_point[0] * loc_point[1]
        sfuncs[4] = -2 * loc_point[1] * (0.5 - loc_point[1])
        sfuncs[5] = 4 * loc_point[1] * (1 - loc_point[0] - loc_point[1])

        return sfuncs

    def eval_shape_gradients(self, loc_point):
        # Evaluate the shape gradients in the local coordinate system
        # Note that no weights are applied!
        sgrads = np.zeros((self._ncount, self._rank))

        sgrads[0, 0] = -3 + 4 * loc_point[0] + 4 * loc_point[1]
        sgrads[0, 1] = -3 + 4 * loc_point[0] + 4 * loc_point[1]
        sgrads[1, 0] = 4 - 8 * loc_point[0] - 4 * loc_point[1]
        sgrads[1, 1] = -4 * loc_point[0]
        sgrads[2, 0] = -1 + 4 * loc_point[0]
        sgrads[2, 1] = 0.0
        sgrads[3, 0] = 4 * loc_point[1]
        sgrads[3, 1] = 4 * loc_point[0]
        sgrads[4, 0] = 0.0
        sgrads[4, 1] = -1 + 4 * loc_point[1]
        sgrads[5, 0] = -4 * loc_point[1]
        sgrads[5, 1] = 4 - 4 * loc_point[0] - 8 * loc_point[1]

        return sgrads


class Line2Shape(Shape):
    def __init__(self, intscheme):
        print('Creating Line2Shape...')

        # Set the nodecount and rank of the elements
        self._ncount = 2
        self._rank = 1

        # Refer to the Shape class to handle the rest of the initialization
        super().__init__(intscheme)

    def get_local_node_coords(self):
        # Return the standard line with nodes at (-1) and (1)
        loc_coords = np.zeros((self._rank, self._ncount))

        loc_coords[0, 0] = -1.0
        loc_coords[0, 1] = 1.0

        return loc_coords

    def eval_shape_functions(self, loc_point):
        # Evalulate the shape functions in the local coordinate system
        sfuncs = np.zeros(self._ncount)

        sfuncs[0] = 0.5 - 0.5 * loc_point[0]
        sfuncs[1] = 0.5 + 0.5 * loc_point[0]

        return sfuncs

    def eval_shape_gradients(self, loc_point):
        # Evaluate the shape gradients in the local coordinate system
        # Note that no weights are applied!
        sgrads = np.zeros((self._ncount, self._rank))

        sgrads[0, 0] = -0.5
        sgrads[1, 0] = 0.5

        return sgrads


class Line3Shape(Shape):
    def __init__(self, intscheme):
        print('Creating Line2Shape...')

        # Set the nodecount and rank of the elements
        self._ncount = 3
        self._rank = 1

        # Refer to the Shape class to handle the rest of the initialization
        super().__init__(intscheme)

    def get_local_node_coords(self):
        # Return the standard line with nodes at (-1) and (1)
        loc_coords = np.zeros((self._rank, self._ncount))

        loc_coords[0, 0] = -1.0
        loc_coords[0, 1] = 0.0
        loc_coords[0, 2] = 1.0

        return loc_coords

    def eval_shape_functions(self, loc_point):
        # Evalulate the shape functions in the local coordinate system
        sfuncs = np.zeros(self._ncount)

        sfuncs[0] = 0.5 * loc_point[0] * (loc_point[0] - 1)
        sfuncs[1] = 1 - loc_point[0]**2
        sfuncs[2] = 0.5 * loc_point[0] * (loc_point[0] + 1)

        return sfuncs

    def eval_shape_gradients(self, loc_point):
        # Evaluate the shape gradients in the local coordinate system
        # Note that no weights are applied!
        sgrads = np.zeros((self._ncount, self._rank))

        sgrads[0, 0] = loc_point[0] - 0.5
        sgrads[1, 0] = -2 * loc_point[0]
        sgrads[2, 0] = loc_point[0] + 0.5

        return sgrads


def declare(factory):
    factory.declare_shape('Triangle3', Tri3Shape)
    factory.declare_shape('Triangle6', Tri6Shape)
    factory.declare_shape('Line2', Line2Shape)
    factory.declare_shape('Line3', Line3Shape)

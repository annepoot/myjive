import numpy as np

from .shape import Shape

__all__ = [
    "Tri3Shape",
    "Tri6Shape",
    "Quad4Shape",
    "Quad9Shape",
    "Line2Shape",
    "Line3Shape",
]


class Tri3Shape(Shape):
    def __init__(self, intscheme):
        # Set the nodecount and rank of the elements
        self._ncount = 3
        self._rank = 2

        # Refer to the Shape class to handle the rest of the initialization
        super().__init__(intscheme)

    @classmethod
    def declare(cls, factory):
        factory.declare_shape("Triangle3", cls)

    def get_local_node_coords(self):
        # Return the standard triangle with nodes at (0,0), (1,0) and (0,1)
        loc_coords = np.zeros((self._ncount, self._rank))

        loc_coords[0, 0] = 0.0
        loc_coords[1, 0] = 1.0
        loc_coords[2, 0] = 0.0
        loc_coords[0, 1] = 0.0
        loc_coords[1, 1] = 0.0
        loc_coords[2, 1] = 1.0

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
        sgrads = np.zeros((self._rank, self._ncount))

        sgrads[0, 0] = -1.0
        sgrads[1, 0] = -1.0
        sgrads[0, 1] = 1.0
        sgrads[1, 1] = 0.0
        sgrads[0, 2] = 0.0
        sgrads[1, 2] = 1.0

        return sgrads

    def get_local_point(self, glob_point, glob_coords):
        # Return the local coordinates corresponding to the given global point
        Ax = glob_coords[0, 0]
        Ay = glob_coords[0, 1]
        Bx = glob_coords[1, 0]
        By = glob_coords[1, 1]
        Cx = glob_coords[2, 0]
        Cy = glob_coords[2, 1]

        mat = np.zeros((self._rank, self._rank))
        rhs = np.zeros(self._rank)

        mat[0, 0] = Bx - Ax
        mat[0, 1] = Cx - Ax
        mat[1, 0] = By - Ay
        mat[1, 1] = Cy - Ay

        rhs[0] = glob_point[0] - Ax
        rhs[1] = glob_point[1] - Ay

        return np.linalg.solve(mat, rhs)

    def contains_local_point(self, loc_point, tol=0.0):
        # Return whether or not the local point falls inside or on the element boundaries
        if loc_point[0] + loc_point[1] > 1 + tol:
            return False
        elif loc_point[0] < 0 - tol:
            return False
        elif loc_point[1] < 0 - tol:
            return False
        else:
            return True

    def contains_global_point(self, glob_point, glob_coords, tol=0.0):
        # Return whether or not the global point falls inside or on the element boundaries
        for i in range(self._ncount):
            ax, ay = glob_coords[i % self._ncount]
            bx, by = glob_coords[(i + 1) % self._ncount]
            px, py = glob_point

            cross = (px - bx) * (ay - by) - (py - by) * (ax - bx)

            if cross <= tol:
                return False
        else:
            return True


class Tri6Shape(Shape):
    def __init__(self, intscheme):
        # Set the nodecount and rank of the elements
        self._ncount = 6
        self._rank = 2

        # Refer to the Shape class to handle the rest of the initialization
        super().__init__(intscheme)

    @classmethod
    def declare(cls, factory):
        factory.declare_shape("Triangle6", cls)

    def get_local_node_coords(self):
        # Return the standard triangle with nodes at (0,0), (1,0) and (0,1)
        loc_coords = np.zeros((self._ncount, self._rank))

        loc_coords[0, 0] = 0.0
        loc_coords[1, 0] = 1.0
        loc_coords[2, 0] = 0.0
        loc_coords[3, 0] = 0.5
        loc_coords[4, 0] = 0.5
        loc_coords[5, 0] = 0.0
        loc_coords[0, 1] = 0.0
        loc_coords[1, 1] = 0.0
        loc_coords[2, 1] = 1.0
        loc_coords[3, 1] = 0.0
        loc_coords[4, 1] = 0.5
        loc_coords[5, 1] = 0.5

        return loc_coords

    def eval_shape_functions(self, loc_point):
        # Evalulate the shape functions in the local coordinate system
        sfuncs = np.zeros(self._ncount)

        sfuncs[0] = (
            2 * (0.5 - loc_point[0] - loc_point[1]) * (1 - loc_point[0] - loc_point[1])
        )
        sfuncs[1] = -2 * loc_point[0] * (0.5 - loc_point[0])
        sfuncs[2] = -2 * loc_point[1] * (0.5 - loc_point[1])
        sfuncs[3] = 4 * loc_point[0] * (1 - loc_point[0] - loc_point[1])
        sfuncs[4] = 4 * loc_point[0] * loc_point[1]
        sfuncs[5] = 4 * loc_point[1] * (1 - loc_point[0] - loc_point[1])

        return sfuncs

    def eval_shape_gradients(self, loc_point):
        # Evaluate the shape gradients in the local coordinate system
        # Note that no weights are applied!
        sgrads = np.zeros((self._rank, self._ncount))

        sgrads[0, 0] = -3 + 4 * loc_point[0] + 4 * loc_point[1]
        sgrads[1, 0] = -3 + 4 * loc_point[0] + 4 * loc_point[1]
        sgrads[0, 1] = -1 + 4 * loc_point[0]
        sgrads[1, 1] = 0.0
        sgrads[0, 2] = 0.0
        sgrads[1, 2] = -1 + 4 * loc_point[1]
        sgrads[0, 3] = 4 - 8 * loc_point[0] - 4 * loc_point[1]
        sgrads[1, 3] = -4 * loc_point[0]
        sgrads[0, 4] = 4 * loc_point[1]
        sgrads[1, 4] = 4 * loc_point[0]
        sgrads[0, 5] = -4 * loc_point[1]
        sgrads[1, 5] = 4 - 4 * loc_point[0] - 8 * loc_point[1]

        return sgrads

    def contains_local_point(self, loc_point, tol=0.0):
        # Return whether or not the local point falls inside or on the element boundaries
        if loc_point[0] + loc_point[1] > 1 + tol:
            return False
        elif loc_point[0] < 0 - tol:
            return False
        elif loc_point[1] < 0 - tol:
            return False
        else:
            return True


class Quad4Shape(Shape):
    def __init__(self, intscheme):
        # Set the nodecount and rank of the elements
        self._ncount = 4
        self._rank = 2

        # Refer to the Shape class to handle the rest of the initialization
        super().__init__(intscheme)

    def get_local_node_coords(self):
        # Return the standard triangle with nodes at (0,0), (1,0) and (0,1)
        loc_coords = np.zeros((self._ncount, self._rank))

        loc_coords[0, 0] = -1.0
        loc_coords[1, 0] = 1.0
        loc_coords[2, 0] = 1.0
        loc_coords[3, 0] = -1.0
        loc_coords[0, 1] = -1.0
        loc_coords[1, 1] = -1.0
        loc_coords[2, 1] = 1.0
        loc_coords[3, 1] = 1.0

        return loc_coords

    def eval_shape_functions(self, loc_point):
        # Evalulate the shape functions in the local coordinate system
        sfuncs = np.zeros(self._ncount)

        x = loc_point[0]
        y = loc_point[1]

        sfuncs[0] = 0.25 * (1 - x) * (1 - y)
        sfuncs[1] = 0.25 * (1 + x) * (1 - y)
        sfuncs[2] = 0.25 * (1 + x) * (1 + y)
        sfuncs[3] = 0.25 * (1 - x) * (1 + y)

        return sfuncs

    def eval_shape_gradients(self, loc_point):
        # Evaluate the shape gradients in the local coordinate system
        # Note that no weights are applied!
        sgrads = np.zeros((self._rank, self._ncount))

        x = loc_point[0]
        y = loc_point[1]

        sgrads[0, 0] = -0.25 * (1 - y)
        sgrads[1, 0] = -0.25 * (1 - x)
        sgrads[0, 1] = 0.25 * (1 - y)
        sgrads[1, 1] = -0.25 * (1 + x)
        sgrads[0, 2] = 0.25 * (1 + y)
        sgrads[1, 2] = 0.25 * (1 + x)
        sgrads[0, 3] = -0.25 * (1 + y)
        sgrads[1, 3] = 0.25 * (1 - x)

        return sgrads

    def contains_local_point(self, loc_point, tol=0.0):
        # Return whether or not the local point falls inside or on the element boundaries
        if loc_point[0] < -1 - tol:
            return False
        elif loc_point[0] > 1 + tol:
            return False
        elif loc_point[1] < -1 - tol:
            return False
        elif loc_point[1] > 1 + tol:
            return False
        else:
            return True

    def contains_global_point(self, glob_point, glob_coords, tol=0.0):
        # Return whether or not the global point falls inside or on the element boundaries
        for i in range(self._ncount):
            ax, ay = glob_coords[i % self._ncount]
            bx, by = glob_coords[(i + 1) % self._ncount]
            px, py = glob_point

            cross = (px - bx) * (ay - by) - (py - by) * (ax - bx)

            if cross < tol:
                return False
        else:
            return True


class Quad9Shape(Shape):
    def __init__(self, intscheme):
        # Set the nodecount and rank of the elements
        self._ncount = 9
        self._rank = 2

        # Refer to the Shape class to handle the rest of the initialization
        super().__init__(intscheme)

    def get_local_node_coords(self):
        # Return the standard triangle with nodes at (0,0), (1,0) and (0,1)
        loc_coords = np.zeros((self._ncount, self._rank))

        loc_coords[0, 0] = -1.0
        loc_coords[1, 0] = 1.0
        loc_coords[2, 0] = 1.0
        loc_coords[3, 0] = -1.0
        loc_coords[4, 0] = 0.0
        loc_coords[5, 0] = 1.0
        loc_coords[6, 0] = 0.0
        loc_coords[7, 0] = -1.0
        loc_coords[8, 0] = 0.0
        loc_coords[0, 1] = -1.0
        loc_coords[1, 1] = -1.0
        loc_coords[2, 1] = 1.0
        loc_coords[3, 1] = 1.0
        loc_coords[4, 1] = -1.0
        loc_coords[5, 1] = 0.0
        loc_coords[6, 1] = 1.0
        loc_coords[7, 1] = 0.0
        loc_coords[8, 1] = 0.0

        return loc_coords

    def eval_shape_functions(self, loc_point):
        # Evalulate the shape functions in the local coordinate system
        sfuncs = np.zeros(self._ncount)

        x = loc_point[0]
        y = loc_point[1]

        sfuncs[0] = 0.25 * x * (x - 1) * y * (y - 1)
        sfuncs[1] = 0.25 * x * (x + 1) * y * (y - 1)
        sfuncs[2] = 0.25 * x * (x + 1) * y * (y + 1)
        sfuncs[3] = 0.25 * x * (x - 1) * y * (y + 1)

        sfuncs[4] = 0.5 * (1 - x**2) * y * (y - 1)
        sfuncs[5] = 0.5 * x * (x + 1) * (1 - y**2)
        sfuncs[6] = 0.5 * (1 - x**2) * y * (y + 1)
        sfuncs[7] = 0.5 * x * (x - 1) * (1 - y**2)

        sfuncs[8] = (1 - x**2) * (1 - y**2)

        return sfuncs

    def eval_shape_gradients(self, loc_point):
        # Evaluate the shape gradients in the local coordinate system
        # Note that no weights are applied!
        sgrads = np.zeros((self._rank, self._ncount))

        x = loc_point[0]
        y = loc_point[1]

        sgrads[0, 0] = 0.5 * (x - 0.5) * y * (y - 1)
        sgrads[1, 0] = 0.5 * x * (x - 1) * (y - 0.5)
        sgrads[0, 1] = 0.5 * (x + 0.5) * y * (y - 1)
        sgrads[1, 1] = 0.5 * x * (x + 1) * (y - 0.5)
        sgrads[0, 2] = 0.5 * (x + 0.5) * y * (y + 1)
        sgrads[1, 2] = 0.5 * x * (x + 1) * (y + 0.5)
        sgrads[0, 3] = 0.5 * (x - 0.5) * y * (y + 1)
        sgrads[1, 3] = 0.5 * x * (x - 1) * (y + 0.5)

        sgrads[0, 4] = -x * y * (y - 1)
        sgrads[1, 4] = (1 - x**2) * (y - 0.5)
        sgrads[0, 5] = (x + 0.5) * (1 - y**2)
        sgrads[1, 5] = x * (x + 1) * -y
        sgrads[0, 6] = -x * y * (y + 1)
        sgrads[1, 6] = (1 - x**2) * (y + 0.5)
        sgrads[0, 7] = (x - 0.5) * (1 - y**2)
        sgrads[1, 7] = x * (x - 1) * -y

        sgrads[0, 8] = (-2 * x) * (1 - y**2)
        sgrads[1, 8] = (1 - x**2) * (-2 * y)

        return sgrads

    def contains_local_point(self, loc_point, tol=0.0):
        # Return whether or not the local point falls inside or on the element boundaries
        if loc_point[0] < -1 - tol:
            return False
        elif loc_point[0] > 1 + tol:
            return False
        elif loc_point[1] < -1 - tol:
            return False
        elif loc_point[1] > 1 + tol:
            return False
        else:
            return True


class Line2Shape(Shape):
    def __init__(self, intscheme):
        # Set the nodecount and rank of the elements
        self._ncount = 2
        self._rank = 1

        # Refer to the Shape class to handle the rest of the initialization
        super().__init__(intscheme)

    def get_local_node_coords(self):
        # Return the standard line with nodes at (-1) and (1)
        loc_coords = np.zeros((self._ncount, self._rank))

        loc_coords[0, 0] = -1.0
        loc_coords[1, 0] = 1.0

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
        sgrads = np.zeros((self._rank, self._ncount))

        sgrads[0, 0] = -0.5
        sgrads[0, 1] = 0.5

        return sgrads

    def get_local_point(self, glob_point, glob_coords):
        # Return the local coordinates corresponding to the given global point
        loc_point = np.zeros(self._rank)

        A = glob_coords[0, 0]
        B = glob_coords[1, 0]
        X = glob_point[0]

        loc_point[0] = (A + B - 2 * X) / (A - B)

        return loc_point

    def contains_local_point(self, loc_point, tol=0.0):
        # Return whether or not the local point falls inside or on the element boundaries
        if loc_point[0] > 1 + tol:
            return False
        elif loc_point[0] < -1 - tol:
            return False
        else:
            return True

    def contains_global_point(self, glob_point, glob_coords, tol=0.0):
        # Return whether or not the global point falls inside or on the element boundaries
        if glob_coords[0, 0] - tol <= glob_point[0] <= glob_coords[1, 0] + tol:
            return True
        else:
            return False


class Line3Shape(Shape):
    def __init__(self, intscheme):
        # Set the nodecount and rank of the elements
        self._ncount = 3
        self._rank = 1

        # Refer to the Shape class to handle the rest of the initialization
        super().__init__(intscheme)

    def get_local_node_coords(self):
        # Return the standard line with nodes at (-1) and (1)
        loc_coords = np.zeros((self._ncount, self._rank))

        loc_coords[0, 0] = -1.0
        loc_coords[1, 0] = 0.0
        loc_coords[2, 0] = 1.0

        return loc_coords

    def eval_shape_functions(self, loc_point):
        # Evalulate the shape functions in the local coordinate system
        sfuncs = np.zeros(self._ncount)

        sfuncs[0] = 0.5 * loc_point[0] * (loc_point[0] - 1)
        sfuncs[1] = 1 - loc_point[0] ** 2
        sfuncs[2] = 0.5 * loc_point[0] * (loc_point[0] + 1)

        return sfuncs

    def eval_shape_gradients(self, loc_point):
        # Evaluate the shape gradients in the local coordinate system
        # Note that no weights are applied!
        sgrads = np.zeros((self._rank, self._ncount))

        sgrads[0, 0] = loc_point[0] - 0.5
        sgrads[0, 1] = -2 * loc_point[0]
        sgrads[0, 2] = loc_point[0] + 0.5

        return sgrads

    def contains_local_point(self, loc_point, tol=0.0):
        # Return whether or not the local point falls inside or on the element boundaries
        if loc_point[0] > 1 + tol:
            return False
        elif loc_point[0] < -1 - tol:
            return False
        else:
            return True

    def contains_global_point(self, glob_point, glob_coords, tol=0.0):
        # Return whether or not the global point falls inside or on the element boundaries
        if glob_coords[0, 0] - tol <= glob_point[0] <= glob_coords[2, 0] + tol:
            return True
        else:
            return False

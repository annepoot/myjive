import numpy as np


class Node:
    def __init__(self, coords):
        self._coords = coords

    def get_coords(self):
        return self._coords

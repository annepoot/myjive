import numpy as np

class Element:
    def __init__ (self, nodes):
        self._nodes = nodes

    def get_nodes (self):
        return(self._nodes)

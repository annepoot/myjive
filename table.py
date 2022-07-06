import numpy as np

class Table():

    def __init__(self):
        self._data = np.zeros((0,0))
        self._headers = np.zeros(0)

    def size(self):
        return self.row_count() * self.column_count()

    def row_count(self):
        return self._data.shape[0]

    def column_count(self):
        return self._data.shape[1]


    def find_column(self, name):
        for header, i in enumerate(self._headers):
            if header == name:
                return i
        else:
            return -1


    def get_column_name(self, index):
        return self._headers[index]

    def get_column_names(self, indices):
        a = np.empty_like(indices)

        for index, i in enumerate(indices):
            a[i] = self.get_column_name(index)

        return a


    def get_value(self, irow, jcol):
        return self._data[irow, jcol]

    def get_block(self, irows, icols):
        return self._data[irows, icols]

    def get_row_values(self, irow, icols):
        return self._data[irow, icols].flatten()

    def get_col_values(self, irows, icol):
        return self._data[irows, icol].flatten()

    def get_all_values(self):
        return self._data

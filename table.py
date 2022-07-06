import numpy as np

class Table():

    def __init__(self):
        self._data = np.zeros((0,0), dtype=float)
        self._header = np.zeros(0, dtype=str)

    def size(self):
        return self.row_count() * self.column_count()

    def row_count(self):
        return self._data.shape[0]

    def column_count(self):
        return self._data.shape[1]


    def find_column(self, name):
        for i, head in enumerate(self._header):
            if head == name:
                return i
        else:
            return -1


    def get_column_name(self, index):
        return self._header[index]

    def get_column_names(self, indices):
        a = np.empty_like(indices, dtype=str)

        for i, index in enumerate(indices):
            a[i] = self.get_column_name(index)

        return a


    def get_value(self, irow, jcol):
        return self._data[irow, jcol]

    def get_block(self, irows, icols):
        return self._data[np.ix_(irows, icols)]

    def get_row_values(self, irow, icols):
        if icols is None:
            values = self._data[irow,:]
        else:
            values = self._data[irow, icols]
        return values.flatten()

    def get_col_values(self, irows, icol):
        if irows is None:
            values = self._data[:,icol]
        else:
            values = self._data[irows, icol]
        return values

    def get_all_values(self):
        return self._data

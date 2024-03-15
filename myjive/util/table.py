import numpy as np

__all__ = ["Table"]


class Table:
    def __init__(self, tbl=None, size=None):
        if tbl is None:
            if size is None:
                self._data = np.zeros((0, 0), dtype=float)
            else:
                self._data = np.zeros((size, 0), dtype=float)
            self._header = np.zeros(0, dtype=str)
        else:
            self._data = tbl._data
            self._header = tbl._header

    def __contains__(self, item):
        return item in self._header

    def __getitem__(self, name):
        jcol = self.get_column(name)
        return self.get_col_values(None, jcol)

    def __iter__(self):
        return iter(self._header)

    def __next__(self):
        return next(self._header)

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

    def find_columns(self, names):
        a = np.empty_like(names, dtype=int)
        for i, name in enumerate(names):
            a[i] = self.find_column(name)
        return a

    def get_column(self, name):
        loc = self.find_column(name)
        if loc < 0:
            raise KeyError("{} could not be found in the table headers".format(name))
        return loc

    def get_columns(self, names):
        a = np.empty_like(names, dtype=int)
        for i, name in enumerate(names):
            a[i] = self.get_column(name)
        return a

    def get_column_name(self, index):
        return self._header[index]

    def get_column_names(self, indices=None):
        if indices is None:
            indices = np.arange(self.column_count())
        a = []
        for index in indices:
            a.append(self.get_column_name(index))
        a = np.array(a)
        return a

    def get_value(self, irow, jcol):
        return self._data[irow, jcol]

    def get_block(self, irows, jcols):
        return self._data[np.ix_(irows, jcols)]

    def get_row_values(self, irow, jcols):
        if jcols is None:
            values = self._data[irow, :]
        else:
            values = self._data[irow, jcols]
        return values.flatten()

    def get_col_values(self, irows, jcol):
        if irows is None:
            values = self._data[:, jcol]
        else:
            values = self._data[irows, jcol]
        return values

    def get_all_values(self):
        return self._data

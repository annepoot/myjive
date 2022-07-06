import numpy as np

from table import Table

class XTable(Table):

    def clear_data(self):
        self._data = np.zeros((0, self._headers.size()))

    def clear_all(self):
        self._header = np.zeros(0)
        self.clear_data()


    def reserve(self, rowcount):
        if self.row_count() < rowcount:
            self._data.resize((rowcount, self.column_count))


    def add_column(self, name):
        if self.find_column(name) < 0:
            self._header = np.append(self._header, name)
            self._data.resize((self.row_count(), self.column_count() + 1))

    def add_columns(self, names):
        for name in names:
            self.add_column(name)


    def set_value(self, irow, jcol, value):
        self._data[irow, jcol] = value

    def add_value(self, irow, jcol, value):
        self.reserve(irow)
        self.set_value(irow, jcol, value)

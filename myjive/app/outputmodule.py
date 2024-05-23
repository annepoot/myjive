import os
import numpy as np
import pandas as pd

from .module import Module
from ..names import GlobNames as gn
from ..util import Table
from ..util.proputils import check_list, get_recursive, split_key

__all__ = ["OutputModule"]


class OutputModule(Module):
    @Module.save_config
    def configure(
        self, globdat, *, files=[r"disp{t}.csv"], keys=[r"state0"], overwrite=False
    ):
        # Validate input arguments
        check_list(self, files)
        check_list(self, keys)
        self._files = files
        self._keys = keys

        self._overwrite = overwrite

        if len(files) != len(keys):
            raise ValueError(
                "'files' and 'values' must have the same number of elements"
            )

    def init(self, globdat):
        pass

    def run(self, globdat):
        for file, key in zip(self._files, self._keys):
            if isinstance(key, list):
                for k in key:
                    header = k.removeprefix("tables.")
                    if "." in k:
                        value = get_recursive(globdat, split_key(k))
                    else:
                        value = globdat[k]
                    lst = self._recursive_output([], value, header)
            else:
                header = key.removeprefix("tables.")
                if "." in key:
                    value = get_recursive(globdat, split_key(key))
                else:
                    value = globdat[key]
                lst = self._recursive_output([], value, header)

            fname = file.format(t=globdat[gn.TIMESTEP])
            if os.path.isfile(fname):
                if self._overwrite:
                    os.remove(fname)
                else:
                    raise RuntimeError("'{}' already exists!".format(fname))

            self._write_to_file(lst, fname)

        return "ok"

    def shutdown(self, globdat):
        pass

    def _recursive_output(self, lst, value, header):
        if isinstance(value, dict):
            for key, val in value.items():
                new_header = self._extend_header(header, key)
                lst = self._recursive_output(lst, val, new_header)
        elif isinstance(value, Table):
            for key in value.get_column_names():
                val = value[key]
                new_header = self._extend_header(header, key)
                lst = self._recursive_output(lst, val, new_header)
        elif isinstance(value, list):
            if hasattr(value, "__len__"):
                for i, val in enumerate(value, 1):
                    new_header = self._extend_header(header, i)
                    lst = self._recursive_output(lst, val, new_header)
            else:
                lst = self._append_single_column(lst, value, header)

        elif isinstance(value, np.ndarray):
            ndim = len(value.shape)
            if ndim == 1:
                lst = self._append_single_column(lst, value, header)
            elif ndim == 2:
                for i, val in enumerate(value):
                    new_header = self._extend_header(header, i)
                    lst = self._recursive_output(lst, val, new_header)
            else:
                raise ValueError("Cannot handle >2D arrays")
        else:
            raise ValueError("Unknown data type")

        return lst

    def _append_single_column(self, lst, value, header):
        lst.append(pd.Series(data=value, name=header))
        return lst

    def _write_to_file(self, lst, fname):
        df = pd.concat(lst, axis=1)

        path = os.path.split(fname)[0]
        if len(path) > 0 and not os.path.isdir(path):
            os.makedirs(path)

        df.to_csv(fname, index=False)

    def _extend_header(self, header, extension):
        new_header = ".".join([header, str(extension)])
        new_header = new_header.replace("..", ".")
        new_header = new_header.removesuffix(".")
        return new_header

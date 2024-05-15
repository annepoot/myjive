import os
import numpy as np

from .module import Module
from ..names import GlobNames as gn
from ..util import Table
from ..util.proputils import check_list, get_recursive

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
        if self._overwrite:
            for file in self._files:
                fname = file.format(t=globdat[gn.TIMESTEP])
                if os.path.isfile(fname):
                    os.remove(fname)

        for file, key in zip(self._files, self._keys):
            fname = file.format(t=globdat[gn.TIMESTEP])

            if isinstance(key, list):
                for k in key:
                    header = k.removeprefix("tables.")
                    if "." in k:
                        value = get_recursive(globdat, self._split_key(k))
                    else:
                        value = globdat[k]
                    self._recursive_output(fname, value, header)
            else:
                header = key.removeprefix("tables.")
                if "." in key:
                    value = get_recursive(globdat, self._split_key(key))
                else:
                    value = globdat[key]
                self._recursive_output(fname, value, header)

        return "ok"

    def shutdown(self, globdat):
        pass

    def _recursive_output(self, fname, value, header):
        if isinstance(value, dict):
            for key, val in value.items():
                new_header = self._extend_header(header, key)
                self._recursive_output(fname, val, new_header)
        elif isinstance(value, Table):
            for key in value.get_column_names():
                val = value[key]
                new_header = self._extend_header(header, key)
                self._recursive_output(fname, val, new_header)
        elif isinstance(value, list):
            if hasattr(value, "__len__"):
                for i, val in enumerate(value, 1):
                    new_header = self._extend_header(header, i)
                    self._recursive_output(fname, val, new_header)
            else:
                self._append_single_column(fname, value, header)

        elif isinstance(value, np.ndarray):
            ndim = len(value.shape)
            if ndim == 1:
                self._append_single_column(fname, value, header)
            elif ndim == 2:
                for i, val in enumerate(value, 1):
                    new_header = self._extend_header(header, i)
                    self._recursive_output(fname, val, new_header)
            else:
                raise ValueError("Cannot handle >2D arrays")
        else:
            raise ValueError("Unknown data type")

    def _append_single_column(self, fname, value, header):
        if os.path.isfile(fname):
            with open(fname, "r") as f:
                lines = f.readlines()

            if len(lines) != len(value) + 1:
                raise ValueError("Incompatible column sizes")

            lines[0] = lines[0].removesuffix("\n") + "," + header + "\n"
            for i, val in enumerate(value, 1):
                lines[i] = lines[i].removesuffix("\n") + "," + str(val) + "\n"

        else:
            lines = [header + "\n"]
            for val in value:
                lines.append(str(val) + "\n")

        path = os.path.split(fname)[0]
        if len(path) > 0 and not os.path.isdir(path):
            os.makedirs(path)

        with open(fname, "w") as f:
            f.writelines(lines)

    def _extend_header(self, header, extension):
        new_header = ".".join([header, str(extension)])
        new_header = new_header.replace("..", ".")
        new_header = new_header.removesuffix(".")
        return new_header

    def _split_key(self, key):
        keys = key.split(".")
        for i, k in enumerate(keys):
            if k.isdigit():
                keys[i] = int(k)
        return keys

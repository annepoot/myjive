from myjive.util import saveconfig as sg

__all__ = ["Material"]


class Material:
    def __init__(self, name):
        self._name = name
        self._config = {}

    def configure(self, globdat):
        raise NotImplementedError("Empty material configure")

    def get_config(self):
        if len(self._config) == 0:
            raise NotImplementedError("Empty material get_config")
        else:
            return self._config

    def get_name(self):
        return self._name

    @classmethod
    def get_type(cls):
        typ = cls.__name__
        if typ[-8:] == "Material":
            typ = typ[:-8]
        return typ

    @staticmethod
    def save_config(configure):
        return sg.save_config(configure)

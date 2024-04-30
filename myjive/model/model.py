from ..util import saveconfig as sg

__all__ = ["Model", "ModelFactory"]


class ModelFactory:
    def __init__(self):
        self._creators = {}

    def declare_model(self, typ, creator):
        self._creators[typ] = creator

    def get_model(self, typ, name):
        creator = self._creators.get(typ)
        if not creator:
            raise ValueError(typ)
        return creator(name)

    def is_model(self, typ):
        return typ in self._creators


class Model:
    def __init__(self, name):
        self._name = name
        self._config = {}

    @classmethod
    def get_type(cls):
        typ = cls.__name__
        if typ[-5:] == "Model":
            typ = typ[:-5]
        return typ

    @classmethod
    def declare(cls, factory):
        typ = cls.get_type()
        factory.declare_model(typ, cls)

    def list_actions(self):
        action_list = []
        for func in dir(self):
            if callable(getattr(self, func)):
                if not func.startswith("_") and func == func.upper():
                    action_list.append(func)
        return action_list

    def take_action(self, action, params, globdat):
        raise (NotImplementedError, "take_action has been deprecated!")

    def configure(self, globdat):
        raise NotImplementedError("Empty model configure")

    def get_config(self):
        if len(self._config) == 0:
            raise NotImplementedError("Empty model get_config")
        else:
            return self._config

    def get_name(self):
        return self._name

    @staticmethod
    def save_config(configure):
        return sg.save_config(configure)

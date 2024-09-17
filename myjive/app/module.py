from ..util import saveconfig as sg

__all__ = ["Module", "ModuleFactory"]


class ModuleFactory:
    def __init__(self):
        self._creators = {}

    def declare_module(self, typ, creator):
        self._creators[typ] = creator

    def get_module(self, typ, name):
        creator = self._creators.get(typ)
        if not creator:
            raise ValueError(typ)
        return creator(name)

    def is_module(self, typ):
        return typ in self._creators


class Module:
    def __init__(self, name):
        self._name = name
        self._config = {}
        self._needs_modelprops = False

    @classmethod
    def get_type(cls):
        typ = cls.__name__
        if typ[-6:] == "Module":
            typ = typ[:-6]
        return typ

    @classmethod
    def declare(cls, factory):
        typ = cls.get_type()
        factory.declare_module(typ, cls)

    def get_relevant_models(self, action, models):
        model_list = []
        for model in models.values():
            if action in model.list_actions():
                model_list.append(model)
        return model_list

    def get_unique_relevant_model(self, action, models):
        if isinstance(action, list):
            model = self.get_unique_relevant_model(action[0], models)
            for act in action[1:]:
                if model is not self.get_unique_relevant_model(act, models):
                    raise RuntimeError(
                        "Multiple relevant models found for '{}' actions".format(action)
                    )
        else:
            model_list = self.get_relevant_models(action, models)
            if len(model_list) < 1:
                raise RuntimeError(
                    "No relevant models found for '{}' action".format(action)
                )
            elif len(model_list) > 1:
                raise RuntimeError(
                    "Multiple relevant models found for '{}' action".format(action)
                )
            else:
                model = model_list[0]
        return model

    def configure(self, globdat):
        raise NotImplementedError("Empty module configure")

    def get_config(self):
        if len(self._config) == 0:
            raise NotImplementedError("Empty module get_config")
        else:
            return self._config

    def init(self, globdat):
        raise NotImplementedError("Empty module init")

    def run(self, globdat):
        raise NotImplementedError("Empty module run")

    def shutdown(self, globdat):
        raise NotImplementedError("Empty module shutdown")

    def get_name(self):
        return self._name

    @staticmethod
    def save_config(configure):
        return sg.save_config(configure)

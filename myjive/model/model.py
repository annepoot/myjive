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

    @classmethod
    def declare(cls, factory):
        name = cls.__name__
        if len(name) > 5 and name[-5:] == "Model":
            name = name[:-5]
        factory.declare_model(name, cls)

    def list_actions(self):
        action_list = []
        for func in dir(self):
            if callable(getattr(self, func)):
                if not func.startswith("_") and func == func.upper():
                    action_list.append(func)
        return action_list

    def take_action(self, action, params, globdat):
        raise (NotImplementedError, "take_action has been deprecated!")

    def configure(self, globdat, **props):
        print("Empty model configure")

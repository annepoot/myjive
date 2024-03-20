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

    @classmethod
    def declare(cls, factory):
        name = cls.__name__
        if len(name) > 6 and name[-6:] == "Module":
            name = name[:-6]
        factory.declare_module(name, cls)

    def get_relevant_models(self, action, models):
        model_list = []
        for model in models:
            if action in model.list_actions():
                model_list.append(model)
        return model_list

    def init(self, globdat, **props):
        print("Empty module init")

    def run(self, globdat):
        print("Empty module run")
        return "exit"

    def shutdown(self, globdat):
        print("Empty module shutdown")

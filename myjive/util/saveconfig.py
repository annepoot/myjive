import inspect


def save_config(configure):
    def wrapper(*args, **kwargs):
        if configure.__name__ != "configure":
            raise RuntimeError("saveconfig can only decorate the configure function.")

        config = {}
        sig = inspect.signature(configure)
        for name in sig.parameters:
            param = sig.parameters[name]
            if param.kind is not param.VAR_KEYWORD:
                if param.default is param.empty:
                    config[name] = None
                else:
                    config[name] = param.default

        for key, arg in zip(config.keys(), args):
            config[key] = arg

        for key, arg in kwargs.items():
            config[key] = arg

        config.pop("globdat")
        self = config.pop("self")
        config["type"] = self.get_type()

        self._config.update(config)

        configure(*args, **kwargs)

    return wrapper

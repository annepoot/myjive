from warnings import warn

from ..declare import declare_all
from ..names import GlobNames as gn
from ..util.proputils import split_off_type

__all__ = ["jive"]


def jive(props, extra_declares=[]):
    # Initialize global database, declare models and modules
    globdat = {}

    # Import the jive models and modules
    declare_all(globdat, extra_declares)

    # Build main Module chain
    print("Configuring module chain...")
    modulefac = globdat[gn.MODULEFACTORY]
    globdat[gn.MODULES] = gather_modules(props, modulefac)
    check_modules(props, globdat[gn.MODULES])

    # Configure modules
    for name, module in globdat[gn.MODULES].items():
        moduleprops = props[name]
        typ, moduleprops = split_off_type(moduleprops)
        module.configure(globdat, **moduleprops)

    # Initialize all modules
    print("Initializing modules")
    for module in globdat[gn.MODULES].values():
        if module._needs_modelprops:
            module.init(globdat, modelprops=props[gn.MODEL])
        else:
            module.init(globdat)

    # Run chain until one of the modules ends the computation
    print("Running chain...")

    for module in globdat[gn.MODULES].values():
        if "exit" in module.run(globdat):
            break

    # Run postprocessing routines
    for module in globdat[gn.MODULES].values():
        module.shutdown(globdat)

    print("End of execution")

    return globdat


def gather_modules(props, module_factory):
    if gn.MODULES in props:
        module_names = props[gn.MODULES]
    else:
        raise ValueError("missing 'modules = [...];' in .pro file")

    chain = {}

    for name in module_names:
        # Get the name of each item in the property file
        if "type" in props[name]:
            typ = props[name]["type"]
        else:
            typ = name.title()
            props[name]["type"] = typ

        # If it refers to a module (and not to a model), add it to the chain
        if module_factory.is_module(typ):
            chain[name] = module_factory.get_module(typ, name)
        else:
            raise ValueError("'{}' is not declared as a module".format(typ))

    return chain


def check_modules(props, modules):
    for module_name in props.keys():
        if module_name not in modules and module_name not in [gn.MODULES, gn.MODEL]:
            warning = "module '{}' defined in props, but not in module list"
            warn(warning.format(module_name))

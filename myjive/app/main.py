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
    print("Initializing module chain...")
    modulefac = globdat[gn.MODULEFACTORY]

    chain = []

    for name in props:
        # Get the name of each item in the property file
        if "type" in props[name]:
            typ = props[name]["type"]
        else:
            typ = name.title()
            props[name]["type"] = typ

        # If it refers to a module (and not to a model), add it to the chain
        if modulefac.is_module(typ):
            chain.append(modulefac.get_module(typ, name))

    # Initialize chain
    for module in chain:
        moduleprops = props[module.get_name()]
        typ, moduleprops = split_off_type(moduleprops)
        module.configure(globdat, **moduleprops)
        if module._needs_modelprops:
            module.init(globdat, modelprops=props[gn.MODELS])
        else:
            module.init(globdat)

    # Run chain until one of the modules ends the computation
    print("Running chain...")

    for module in chain:
        if "exit" in module.run(globdat):
            break

    # Run postprocessing routines
    for module in chain:
        module.shutdown(globdat)

    print("End of execution")

    return globdat

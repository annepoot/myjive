import declare
from names import GlobNames as gn


def jive(props):
    # Initialize global database, declare models and modules
    globdat = {}

    declare.declare_models(globdat)
    declare.declare_modules(globdat)
    declare.declare_shapes(globdat)

    # Build main Module chain
    print('Initializing module chain...')
    modulefac = globdat[gn.MODULEFACTORY]
    allmodules = ['Init', 'Solver', 'Nonlin', 'Arclen', 'VTKout', 'FrameView', 'LinBuck', 'Output']

    # Only select modules defined in the .pro file
    modules = [i for i in allmodules if i.lower() in props]
    chain = []
    for module in modules:
        chain.append(modulefac.get_module(module, module.lower()))

    # Initialize chain
    for module in chain:
        module.init(props, globdat)

    # Run chain until one of the modules ends the computation
    print('Running chain...')

    keep_going = True

    while keep_going:
        for module in chain:
            if 'exit' in module.run(globdat):
                keep_going = False

    # Run postprocessing routines
    for module in chain:
        module.shutdown(globdat)

    print('End of execution')

    return globdat

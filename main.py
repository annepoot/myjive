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

    chain = [modulefac.get_module('Init', 'init'),
             modulefac.get_module('Solver', 'solver'),
             modulefac.get_module('VTKOut', 'vtkout'),
             modulefac.get_module('FrameView', 'frameview')]

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

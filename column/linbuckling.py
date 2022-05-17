import sys
sys.path.append('../')

import numpy as np

import proputils as pu
import main

import declare

from names import GlobNames as gn


props = pu.parse_file('column.pro')

globdat = main.jive(props)

declare.declare_models(globdat)
declare.declare_modules(globdat)
declare.declare_shapes(globdat)

print('Initializing module chain...')

modulefac = globdat[gn.MODULEFACTORY]

chain = [modulefac.get_module('Init', 'init'),
         modulefac.get_module('LinBuck', 'linbuck'),
         modulefac.get_module('VTKOut', 'vtkout'),
         modulefac.get_module('FrameView', 'frameview')]

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

import sys

import jive.util.proputils as pu
from jive.app import main

if len(sys.argv) != 2:
    raise RuntimeError('Script expects exactly one argument')

# Read input file

print('Reading model input...')
props = pu.parse_file(sys.argv[1])

# Run Jive

globdat = main.jive(props)

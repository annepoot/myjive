import sys

import fileparser as fp
import main

if len(sys.argv) != 2:
    raise RuntimeError ( 'Script expects exactly one argument' )

# Read input file

print('Reading model input...')
props = fp.parse (sys.argv[1])

# Run Jive

main.jive (props)



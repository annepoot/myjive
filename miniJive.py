import sys

import myjive.util.proputils as pu
from myjive.app import main
from myjivex import declare_all as declarex

if len(sys.argv) != 2:
    raise RuntimeError("Script expects exactly one argument")

# Read input file

print("Reading model input...")
props = pu.parse_file(sys.argv[1])

# Run Jive

globdat = main.jive(props, extra_declares=[declarex])

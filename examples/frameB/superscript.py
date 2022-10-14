import sys
sys.path.append('../../')

import main
import proputils as pu

props = pu.parse_file('frame.pro')
globdat = main.jive(props)

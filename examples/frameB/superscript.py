import sys
sys.path.append('../../')

from jive.app import main
import jive.util.proputils as pu

props = pu.parse_file('frame.pro')
globdat = main.jive(props)

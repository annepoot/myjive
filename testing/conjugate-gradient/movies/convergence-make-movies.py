import sys
sys.path.append('../../../')

import plotutils as pu

preconditioners = ['none', 'diag', 'ichol']

for P in preconditioners:

    modes = ['solutions', 'residuals', 'differences']

    for mode in modes:
        framepath = 'frames/{}/{}'.format(mode, P)
        moviename = '{}-{}.avi'.format(mode, P)

        print('Generating movie: {}'.format(moviename))

        pu.make_movie(framepath, moviename)

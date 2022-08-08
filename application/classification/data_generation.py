import sys
import os

sys.path.append('../')

import numpy as np
import pandas as pd
from datetime import datetime
import main
import proputils as pu
from names import GlobNames as gn

# Read the problem props from the general beam file
props = pu.parse_file('beam.pro')

# Set all general variables beforehand (to also easily get them in the README file)
nsamples = 1
E_mean = np.log(10000)
E_std = 0.1
nu_mean = np.log(0.2)
nu_std = 0.1
n_det_rule = "int(20 * sample / nsamples)"
risk_rule = "'unnecessary' if n_det < 8 else 'maintenance' if n_det < 16 else 'demolition'"
fname = 'classification-data.csv'
locX = props['model']['solid']['material']['locX']
locY = props['model']['solid']['material']['locY']
stdX = props['model']['solid']['material']['stdX']
stdY = props['model']['solid']['material']['stdY']
scale = props['model']['solid']['material']['scale']

# Use a random number generator for the random sampling stuff
rng = np.random.default_rng(0)

# Create an empty dataframe, in which to store all data
rows_list = []

# Go over all samples
for sample in range(nsamples):

    print('\n' + 50 * '=')
    print('\tSAMPLE {} OUT OF {}'.format(sample+1, nsamples))
    print(50 * '=' + '\n')

    # Get E and nu from a random distribution
    E = rng.lognormal(E_mean, E_std)
    nu = rng.lognormal(nu_mean, nu_std)

    # Get the number of deteriorations and intervention from the rules set beforehand
    n_det = eval(n_det_rule, {'sample':sample, 'nsamples':nsamples})
    risk = eval(risk_rule, {'n_det':n_det})

    # Set these variables in the model
    props['model']['solid']['material']['E'] = E
    props['model']['solid']['material']['nu'] = nu
    props['model']['solid']['material']['deteriorations'] = n_det
    props['model']['solid']['material']['seed'] = sample

    globdat = main.jive(props)

    # Move the created vtu file to a separate folder
    cwd = os.getcwd()
    os.replace(cwd + '/stiffness1.vtu', cwd + '/vtk-output/stiffness{}.vtu'.format(sample+1))

    # Get the nodes, displacements and stiffnesses
    nodes = globdat[gn.NSET]
    u = globdat[gn.STATE0]
    dx = u[:len(u)//2]
    dy = u[len(u)//2:]
    E_true = globdat[gn.TABLES]['stiffness']['']

    # Write all relevant outputs to a single file
    for i, node in enumerate(nodes):
        coords = node.get_coords()

        row_dict = {
            'sample': sample,
            'E_pure': E,
            'nu': nu,
            'deteriorations': n_det,
            'intervention': risk,
            'node': i,
            'x': coords[0],
            'y': coords[1],
            'dx': dx[i],
            'dy': dy[i],
            'E_true': E_true[i]
        }

        rows_list.append(row_dict)

# Generate a dataframe and write it to an output csv
df = pd.DataFrame(rows_list)
df.to_csv(fname, index=False)

# Get the time and size of the CSV file that has been generated
time = datetime.fromtimestamp(os.path.getmtime(fname)).strftime('%Y-%m-%d %H:%M:%S')
size = os.path.getsize(fname) / 1024**2

# Write all relevant info to a README file
with open('README.txt', 'w') as f:
    f.write('FILE INFO\n')
    f.write('output file name: {:s}\n'.format(fname))
    f.write('output creation date: {:s}\n'.format(time))
    f.write('output file size: {:.2f} MB\n'.format(size))
    f.write('\nPROBLEM INFO\n')
    f.write('number of samples: {}\n'.format(nsamples))
    f.write('E ~ lognormal(mean: ln({:.0f}), std: {:.1f})\n'.format(np.exp(E_mean), E_std))
    f.write('nu ~ lognormal(mean: ln({:.1f}), std: {:.1f})\n'.format(np.exp(nu_mean), nu_std))
    f.write('n_det = {:s}\n'.format(n_det_rule))
    f.write('intervention = {:s}\n'.format(risk_rule))
    f.write('deterioration note: in the two expressions below, x and y are the center coordinates of a random element\n')
    f.write('deterioration mean: ({}, {})\n'.format(locX, locY))
    f.write('deterioration std: ({}, {})\n'.format(stdX, stdY))
    f.write('deterioration scale: {}\n'.format(scale))

FILE INFO
output file name: classification-data.csv
output creation date: 2022-08-08 17:21:49
output file size: 0.10 MB

PROBLEM INFO
number of samples: 1
E ~ lognormal(mean: ln(10000), std: 0.1)
nu ~ lognormal(mean: ln(0.2), std: 0.1)
n_det = int(20 * sample / nsamples)
intervention = 'unnecessary' if n_det < 8 else 'maintenance' if n_det < 16 else 'demolition'
deterioration note: in the two expressions below, x and y are the center coordinates of a random element
deterioration mean: (x, y)
deterioration std: (np.random.uniform(0.5,1.0), np.random.uniform(0.5,1.0))

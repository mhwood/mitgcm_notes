

import os
import numpy as np


config_dir = '../'

n_rows = 31
n_cols = 90

max_depth = -500

bathy_grid = max_depth*np.ones((n_rows, n_cols))

bathy_grid[:,0] = 0
bathy_grid[:,-1] = 0
bathy_grid[0,:] = 0
bathy_grid[-1,:] = 0

output_file = os.path.join(config_dir,'input','bathymetry.bin')
bathy_grid.ravel(order='C').astype('>f4').tofile(output_file)

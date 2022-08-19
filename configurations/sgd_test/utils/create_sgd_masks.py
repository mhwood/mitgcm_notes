

import os
import numpy as np


config_dir = '../'

n_rows = 31
n_cols = 90
Nr = 50


mask_grid = np.zeros((Nr, n_rows, n_cols))

mask_grid[Nr-2,16,1] = 3
mask_grid[Nr-3,16,1] = 2
mask_grid[Nr-4,16,1] = 1

if 'sgd' not in os.listdir(os.path.join(config_dir,'input')):
    os.mkdir(os.path.join(config_dir,'input','sgd'))

output_file = os.path.join(config_dir,'input','sgd','sgd_outlet_1_mask.bin')
mask_grid.ravel(order='C').astype('>f4').tofile(output_file)

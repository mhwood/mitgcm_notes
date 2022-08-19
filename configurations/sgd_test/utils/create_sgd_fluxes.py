

import os
import numpy as np
import matplotlib.pyplot as plt

config_dir = '../'

n_timesteps = 720
n_points = 3

A49 = 500
A48 = 1000
A47 = 250

flux_grid = np.zeros((n_timesteps,n_points))

flux_grid[:,0] = A47*np.cos(4*np.pi*np.arange(n_timesteps)/n_timesteps) + A47
flux_grid[:,1] = A48*np.cos(4*np.pi*np.arange(n_timesteps)/n_timesteps) + A48
flux_grid[:,2] = A49*np.cos(4*np.pi*np.arange(n_timesteps)/n_timesteps) + A49

plt.plot(flux_grid[:,0])
plt.plot(flux_grid[:,1])
plt.plot(flux_grid[:,2])
plt.show()

if 'sgd' not in os.listdir(os.path.join(config_dir,'input')):
    os.mkdir(os.path.join(config_dir,'input','sgd'))

output_file = os.path.join(config_dir,'input','sgd','sgd_outlet_1_flux.bin')
flux_grid.ravel(order='C').astype('>f4').tofile(output_file)

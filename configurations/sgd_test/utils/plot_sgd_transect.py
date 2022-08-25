

import os
import numpy as np
import netCDF4 as nc4
import matplotlib.pyplot as plt
import cmocean.cm as cm

config_dir = '../'

temp_file = os.path.join(config_dir,'run','mnc_0001','theta.0000000000.t001.nc')
ds = nc4.Dataset(temp_file)
theta = ds.variables['THETA'][:,:,:,:]
ds.close()

salt_file = os.path.join(config_dir,'run','mnc_0001','salt.0000000000.t001.nc')
ds = nc4.Dataset(salt_file)
salt = ds.variables['SALT'][:,:,:,:]
ds.close()

x = np.arange(90)
y = np.arange(31)
z = -1*np.arange(50)

X, Z = np.meshgrid(x,z)

for i in range(72):

    plt.subplot(3,2,1)
    C = plt.pcolormesh(X, Z, theta[i,:,16,:], vmin=0,vmax=1, cmap=cm.thermal)
    plt.colorbar(C)
    plt.title('Temperature')

    plt.subplot(3, 2, 2)
    C = plt.pcolormesh(X, Z, salt[i, :, 16, :], vmin=34.5, vmax=35, cmap = cm.haline)
    plt.colorbar(C)
    plt.title('Salinity')

    plt.subplot(3, 2, 3)
    C = plt.pcolormesh(X, Z, theta[i, :, 16, :] - theta[0, :, 16, :],
                       vmin=-0.1, vmax=0.1,cmap='seismic')
    plt.colorbar(C)
    plt.title('Temperature Anomaly')

    plt.subplot(3, 2, 4)
    C = plt.pcolormesh(X, Z, salt[i, :, 16, :] - salt[0, :, 16, :],
                       vmin=-0.05, vmax=0.05,cmap='seismic')
    plt.colorbar(C)
    plt.title('Salinity Anomaly')

    plt.suptitle('Timestep = '+str(i))

    plt.show()


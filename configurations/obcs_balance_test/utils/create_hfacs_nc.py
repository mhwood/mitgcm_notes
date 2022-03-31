
import os
import simplegrid as sg
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
from pyproj import Transformer
import argparse
import ast
import sys

########################################################################################################################

def create_hfacs_file(config_dir,bathy_file):

    print('Outputting an nc file with the hfac grids')

    n_rows = 20
    n_cols = 20

    bathy_file = os.path.join(config_dir,'input', bathy_file)
    bathy_grid = np.fromfile(bathy_file,'>f4')
    bathy_grid = np.reshape(bathy_grid, (n_rows, n_cols))

    delR = np.array([1.00,    1.14,    1.30,    1.49,   1.70,
                      1.93,    2.20,    2.50,    2.84,   3.21,
                      3.63,    4.10,    4.61,    5.18,   5.79,
                      6.47,    7.20,    7.98,    8.83,   9.73,
                     10.69,   11.70,   12.76,   13.87,  15.03,
                     16.22,   17.45,   18.70,   19.97,  21.27,
                     22.56,   23.87,   25.17,   26.46,  27.74,
                     29.00,   30.24,   31.45,   32.65,  33.82,
                     34.97,   36.09,   37.20,   38.29,  39.37,
                     40.45,   41.53,   42.62,   43.73,  44.87,
                     46.05,   47.28,   48.56,   49.93,  51.38,
                     52.93,   54.61,   56.42,   58.38,  60.53,
                     62.87,   65.43,   68.24,   71.33,  74.73,
                     78.47,   82.61,   87.17,   92.21,  97.79,
                    103.96,  110.79,  118.35,  126.73, 136.01,
                    146.30,  157.71,  170.35,  184.37, 199.89,
                    217.09,  236.13,  257.21,  280.50, 306.24,
                    334.64,  365.93,  400.38,  438.23, 479.74,])

    Z_bottom = np.cumsum(delR)
    Z_top = np.concatenate([np.array([0]),Z_bottom[:-1]])
    Z = (Z_bottom + Z_top)/2

    # sys.path.insert(1, os.path.join(config_dir, 'utils'))
    import downscale_functions as df

    print('   - Calculating the hFacC grid')
    hFac_C = df.create_hFacC_grid(bathy_grid,delR, hFacMin=0.3, hFacMinDr=1.0)
    print('   - Calculating the hFacS grid')
    hFac_S = df.create_hFacS_grid(bathy_grid,delR, hFacMin=0.3, hFacMinDr=1.0)
    print('   - Calculating the hFacW grid')
    hFac_W = df.create_hFacW_grid(bathy_grid,delR, hFacMin=0.3, hFacMinDr=1.0)

    mitgrid_file = os.path.join(config_dir,'input', 'tile001.mitgrid')
    grid_dict = sg.gridio.read_mitgridfile(mitgrid_file, n_cols, n_rows)
    XC = grid_dict['XC'].T
    YC = grid_dict['YC'].T

    output_file = os.path.join('..',  'input', 'hfac_grids.nc')
    ds = nc4.Dataset(output_file,'w')

    ds.createDimension('x',np.shape(XC)[1])
    ds.createDimension('y', np.shape(XC)[0])
    ds.createDimension('z', len(delR))

    ds.source=bathy_file

    var = ds.createVariable('XC','f4',('y','x'))
    var[:] = XC
    var = ds.createVariable('YC', 'f4', ('y','x'))
    var[:] = YC
    var = ds.createVariable('z', 'f4', ('z',))
    var[:] = Z

    var = ds.createVariable('hFacC', 'f4', ('z','y', 'x'))
    var[:,:,:] = hFac_C
    var = ds.createVariable('hFacW', 'f4', ('z','y', 'x'))
    var[:,:,:] = hFac_W
    var = ds.createVariable('hFacS', 'f4', ('z','y', 'x'))
    var[:,:,:] = hFac_S

    ds.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    parser.add_argument("-b", "--bathy_file", action="store",
                        help="The bathymetry file to read in.", dest="bathy_file",
                        type=str, required=False, default='bathymetry.bin')

    args = parser.parse_args()
    config_dir = args.config_dir
    bathy_file = args.bathy_file

    create_hfacs_file(config_dir, bathy_file)

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

def create_hfacs_file(config_dir):

    print('Outputting an nc file with the hfac grids')

    n_rows = 20
    n_cols = 20


    # read in the wet grids made by python
    wetgrid_file = os.path.join('..', 'input', 'hfac_grids.nc')
    ds = nc4.Dataset(wetgrid_file)
    python_hFacC = np.array(ds.variables['hFacC'][:, :, :])
    python_hFacS = np.array(ds.variables['hFacS'][:, :, :])
    python_hFacW = np.array(ds.variables['hFacW'][:, :, :])
    ds.close()

    # read in the wet grids made by the model
    wetgrid_file = os.path.join('..', 'run_unbalanced','mnc_0001', 'grid.t001.nc')
    ds = nc4.Dataset(wetgrid_file)
    mitgcm_hFacC = np.array(ds.variables['HFacC'][:, :, :])
    mitgcm_hFacS = np.array(ds.variables['HFacS'][:, :, :])
    mitgcm_hFacW = np.array(ds.variables['HFacW'][:, :, :])
    ds.close()

    test_row = 10
    test_col = 10
    python_profile = python_hFacC[:,test_row, test_col]
    mitgcm_profile = mitgcm_hFacC[:, test_row, test_col]

    print(np.sum(mitgcm_hFacC!=python_hFacC)),print(np.size(python_hFacC))

    bad_depths, bad_rows, bad_cols = np.where(mitgcm_hFacC!=python_hFacC)
    for i in range(len(bad_depths)):
        python_estimate = python_hFacC[bad_depths[i],bad_rows[i],bad_cols[i]]
        mitgcm_estimate = mitgcm_hFacC[bad_depths[i],bad_rows[i],bad_cols[i]]
        if np.abs(python_estimate-mitgcm_estimate)>0.001:
            print(bad_depths[i],bad_rows[i],bad_cols[i],python_estimate,mitgcm_estimate)

    # plt.plot(python_profile,'b-')
    # plt.plot(mitgcm_profile, 'g-')
    # plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_hfacs_file(config_dir)
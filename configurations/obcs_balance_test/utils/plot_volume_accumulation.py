
import os
import simplegrid as sg
import netCDF4 as nc4
import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import argparse

def read_grid_variables(config_dir):

    wetgrid_file = os.path.join('..', 'input', 'grid.t001.nc')
    ds = nc4.Dataset(wetgrid_file)
    rA = np.array(ds.variables['rA'][:, :])
    ds.close()

    return(rA)

def read_volume_timeseries(config_dir,experiment,rA):

    wetgrid_file = os.path.join('..', 'run_'+experiment, 'mnc_0001', 'ETAN.0000000000.t001.nc')
    ds = nc4.Dataset(wetgrid_file)
    EtaN = np.array(ds.variables['ETAN'][:, :, :, :])
    EtaN = EtaN[:,0,:,:]
    ds.close()

    modeled_volume_timeseries = np.zeros((np.shape(EtaN)[0],))
    modeled_mean_EtaN_timeseries = np.zeros((np.shape(EtaN)[0],))

    for i in range(len(modeled_volume_timeseries)):
        modeled_volume_timeseries[i] = np.sum(EtaN[i, :, :] * rA)
        modeled_mean_EtaN_timeseries[i] = np.mean(EtaN[i, :, :])

    return(modeled_volume_timeseries, modeled_mean_EtaN_timeseries)


def read_expected_timeseries(config_dir):
    ds = nc4.Dataset(os.path.join('..', 'input', 'expected_unbalanced_results.nc'))
    expected_volume_timeseries = ds.variables['volume'][:]
    expected_EtaN_timeseries = ds.variables['mean_etan'][:]
    ds.close()
    return(expected_volume_timeseries, expected_EtaN_timeseries)


########################################################################################################################

def plot_accumulation(config_dir):

    llc = 90

    n_rows = 20
    n_cols = 20

    rA = read_grid_variables(config_dir)

    zero_flux_volume_timeseries, zero_flux_EtaN_timeseries = read_volume_timeseries(config_dir, 'zero_flux',rA)
    unbalanced_volume_timeseries, unbalanced_EtaN_timeseries = read_volume_timeseries(config_dir, 'unbalanced', rA)
    balanced_by_mitgcm_volume_timeseries, balanced_by_mitgcm_EtaN_timeseries = read_volume_timeseries(config_dir, 'balanced_by_mitgcm', rA)
    expected_volume_timeseries, expected_EtaN_timeseries = read_expected_timeseries(config_dir)

    plt.subplot(1, 2, 1)
    plt.plot(zero_flux_volume_timeseries, 'k--', label='zero_flux')
    plt.plot(unbalanced_volume_timeseries, 'b-', label='unbalanced')
    plt.plot(expected_volume_timeseries, 'b--', label='expected imbalance')
    plt.plot(balanced_by_mitgcm_volume_timeseries, 'g-', label='balanced_by_mitgcm')
    plt.legend()
    plt.subplot(1, 2, 2)
    plt.plot(zero_flux_EtaN_timeseries,'k--',label='zero_flux')
    plt.plot(unbalanced_EtaN_timeseries, 'b-', label='unbalanced')
    plt.plot(expected_EtaN_timeseries, 'b--', label='expected imbalance')
    plt.plot(balanced_by_mitgcm_EtaN_timeseries, 'g-', label='balanced_by_mitgcm')
    plt.legend()
    plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    plot_accumulation(config_dir)
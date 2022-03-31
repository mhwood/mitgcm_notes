
import os
import simplegrid as sg
import netCDF4 as nc4
import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import argparse

def read_grid_variables(config_dir):

    wetgrid_file = os.path.join(config_dir, 'input', 'grid.t001.nc')
    ds = nc4.Dataset(wetgrid_file)
    rA = np.array(ds.variables['rA'][:, :])
    ds.close()

    return(rA)

def read_volume_timeseries(config_dir,experiment,rA):

    wetgrid_file = os.path.join(config_dir, 'run_'+experiment, 'mnc_0001', 'ETAN.0000000000.t001.nc')
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

def read_expected_timeseries(config_dir,experiment):
    ds = nc4.Dataset(os.path.join(config_dir, 'input', 'expected_'+experiment+'_results.nc'))
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


    # get the unbalanced timeseries
    unbalanced_constant_volume_timeseries, _ = read_volume_timeseries(config_dir,'unbalanced_constant',rA)
    unbalanced_periodic_volume_timeseries, _ = read_volume_timeseries(config_dir, 'unbalanced_periodic', rA)

    expected_unbalanced_constant_volume_timeseries, _ = read_expected_timeseries(config_dir, 'unbalanced_constant')
    expected_unbalanced_periodic_volume_timeseries, _ = read_expected_timeseries(config_dir,'unbalanced_periodic')

    balanced_constant_volume_timeseries, _ = read_volume_timeseries(config_dir,'balanced_constant_manual',rA)
    balanced_periodic_volume_timeseries, _ = read_volume_timeseries(config_dir,'balanced_periodic_manual',rA)


    plt.subplot(1, 2, 1)
    # plt.plot(expected_unbalanced_constant_volume_timeseries, 'k.', label='expected imbalance', markersize = 10)
    # plt.plot(unbalanced_constant_volume_timeseries, 'b-', label='unbalanced')
    plt.plot(balanced_constant_volume_timeseries, 'g-', label='balanced')
    plt.title('Constant Boundary Imbalance')
    plt.legend()

    plt.subplot(1, 2, 2)
    # plt.plot(expected_unbalanced_periodic_volume_timeseries, 'k.', label='expected imbalance', markersize = 10)
    # plt.plot(unbalanced_periodic_volume_timeseries, 'b-', label='unbalanced')
    plt.plot(balanced_periodic_volume_timeseries, 'g-', label='balanced')
    plt.legend()
    plt.title('Periodic Boundary Imbalance')

    plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L2, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    plot_accumulation(config_dir)

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
    dxG = np.array(ds.variables['dxG'][:, :])
    dyG = np.array(ds.variables['dyG'][:, :])
    hFacS = np.array(ds.variables['HFacS'][:, :, :])
    hFacW = np.array(ds.variables['HFacW'][:, :, :])
    drF = np.array(ds.variables['drF'][:])
    ds.close()

    return(rA, dxG, dyG, drF, hFacS, hFacW)

def get_reference_profiles():
    tRef = np.array([18.89, 18.89, 18.89, 18.89, 18.89, 18.87,
                    18.85, 18.82, 18.80, 18.73, 18.65, 18.57,
                    18.40, 18.22, 18.00, 17.74, 17.44, 17.12,
                    16.76, 16.39, 15.98, 15.55, 15.08, 14.59,
                    14.07, 13.53, 12.99, 12.47, 11.97, 11.49,
                    11.02, 10.57, 10.12, 9.71, 9.27, 8.88,
                    8.46, 8.09, 7.71, 7.37, 7.03, 6.72,
                    6.42, 6.13, 5.86, 5.59, 5.34, 5.09,
                    4.87, 4.65, 4.45, 4.26, 4.08, 3.91,
                    3.75, 3.60, 3.47, 3.33, 3.20, 3.08,
                    2.96, 2.84, 2.73, 2.62, 2.51, 2.42,
                    2.32, 2.23, 2.14, 2.06, 1.98, 1.90,
                    1.81, 1.73, 1.65, 1.57, 1.49, 1.41,
                    1.33, 1.24, 1.15, 1.06, 0.98, 0.94,
                    0.91, 0.92, 0.98, 0.98, 0.98, 0.98])


    sRef = np.array([34.84, 34.84, 34.84, 34.84, 34.84, 34.84,
                    34.85, 34.85, 34.85, 34.86, 34.87, 34.88,
                    34.89, 34.90, 34.92, 34.94, 34.96, 34.98,
                    35.00, 35.02, 35.04, 35.06, 35.07, 35.07,
                    35.07, 35.05, 35.03, 35.01, 34.98, 34.95,
                    34.92, 34.89, 34.85, 34.82, 34.79, 34.76,
                    34.73, 34.71, 34.68, 34.66, 34.64, 34.62,
                    34.61, 34.60, 34.59, 34.59, 34.58, 34.58,
                    34.59, 34.59, 34.60, 34.60, 34.61, 34.62,
                    34.63, 34.64, 34.65, 34.66, 34.67, 34.68,
                    34.69, 34.70, 34.71, 34.71, 34.72, 34.72,
                    34.73, 34.73, 34.74, 34.74, 34.74, 34.74,
                    34.75, 34.74, 34.74, 34.74, 34.74, 34.74,
                    34.74, 34.74, 34.73, 34.73, 34.73, 34.73,
                    34.73, 34.72, 34.72, 34.72, 34.72, 34.72])

    return(tRef,sRef)

def calculate_flux(boundary,field,width,delR,hFac):
    flux = np.zeros_like(field)
    total_area = 0
    for i in range(np.shape(field)[1]):
        for j in range(np.shape(field)[0]):
            flux[j,i] = width[i]*delR[j]*hFac[j,i]*field[j,i]
            total_area += width[i]*delR[j]*hFac[j,i]

            # if i==11 and boundary=='south' and j==82:
            #     print('|- south -----------------------------|')
            #     print('    drF(k) ',delR[j])
            #     print('    hFacS[j,i] ',hFac[j,i])
            #     print('    dxG[i] ',width[i])
            #     print('    OBSv[j,i] ',field[j,i])
            #
            # if i==10 and boundary=='north' and j==68:
            #     print('|- north -----------------------------|')
            #     print('    drF(k) ',delR[j])
            #     print('    hFacS[j,i] ',hFac[j,i])
            #     print('    dxG[i] ',width[i])
            #     print('    OBNv[j,i] ',field[j,i])
            #
            # if i==5 and boundary=='west' and j==85:
            #     print('|- west -----------------------------|')
            #     print('    drF(k) ',delR[j])
            #     print('    hFacW[j,i] ',hFac[j,i])
            #     print('    dyG[i] ',width[i])
            #     print('    OBWu[j,i] ',field[j,i])

    total_flux = np.sum(flux)
    return(total_flux,total_area)

def create_boundary_condition(condition_type, boundary,var_name,
                              dxG, dyG, drF,
                              hFacS, hFacW):
    n_timestep = 8*24
    Nr = 90
    sNx = 20
    sNy = 20

    continue_to_calculation = False
    if boundary=='north':
        sN = sNx
        hFac = hFacS[:,-2,:]
        width = dxG[-2,:]
        if var_name=='VVEL':
            continue_to_calculation = True
    if boundary=='south':
        sN = sNx
        hFac = hFacS[:,1,:]
        width = dxG[1, :]
        if var_name=='VVEL':
            continue_to_calculation = True
    if boundary=='west':
        sN = sNy
        hFac = hFacW[:,:,1]
        width = dyG[:, 1]
        if var_name=='UVEL':
            continue_to_calculation = True

    # note: the first and last grid cells are not calculated in the flux
    # one way to solve this issue is just to make hfac=0 in these areas?
    hFac[:, 0] = 0
    hFac[:, -1] = 0

    grid = np.zeros((n_timestep,Nr,sN))

    flux_timeseries = np.zeros((n_timestep,))

    if condition_type=='unbalanced_constant':

        if continue_to_calculation:
            normalized_vel_profile = 1/drF
            normalized_vel_profile = normalized_vel_profile-normalized_vel_profile[-1]
            normalized_vel_profile = normalized_vel_profile*(1/normalized_vel_profile[0])
            if boundary=='south' or boundary=='north':
                normalized_vel_profile *= 2

            vel_grid = np.zeros((Nr,sN))
            for i in range(sN):
                vel_grid[:,i] = 0.1*normalized_vel_profile
            vel_grid *= hFac

            total_flux, total_area = calculate_flux(boundary,vel_grid,width,drF,hFac)
            flux_timeseries[:] = total_flux

            # print('    area ', total_area)
            # print('    flow ', total_flux)

            for timestep in range(n_timestep):
                grid[timestep,:,:] = vel_grid

    if condition_type=='unbalanced_periodic':

        if continue_to_calculation:
            normalized_vel_profile = 1/drF
            normalized_vel_profile = normalized_vel_profile-normalized_vel_profile[-1]
            normalized_vel_profile = normalized_vel_profile*(1/normalized_vel_profile[0])

            if boundary=='south' or boundary=='north':
                normalized_vel_profile *= 2

            vel_grid = np.zeros((Nr,sN))
            for i in range(sN):
                vel_grid[:,i] = 0.1*normalized_vel_profile
            vel_grid *= hFac

            # print('    area ', total_area)
            # print('    flow ', total_flux)

            amplitude = 2
            period = 12

            for timestep in range(n_timestep):
                grid[timestep,:,:] = vel_grid+ vel_grid*amplitude*np.sin(2*np.pi*timestep/period)
                total_flux, total_area = calculate_flux(boundary, grid[timestep,:,:], width, drF, hFac)
                flux_timeseries[timestep] = total_flux

    return(grid,flux_timeseries)


########################################################################################################################

def create_boundary_conditions(config_dir):

    print('Creating sets of boundary conditions')

    print('    - Reading in the grid information')
    rA, dxG, dyG, drF, hFacS, hFacW = read_grid_variables(config_dir)
    total_area = np.sum(dxG[:-1,:]*dyG[:,:-1])

    for condition_type in ['unbalanced_constant','unbalanced_periodic']:
        print('    - Creating the boundary conditions for the '+condition_type+' experiment')

        if 'obcs_' + condition_type not in os.listdir(os.path.join(config_dir, 'input')):
            os.mkdir(os.path.join(config_dir, 'input', 'obcs_' + condition_type))

        flux_started = False

        for boundary in ['west','north','south']:
            for var_name in ['UVEL','VVEL']:
                bc_grid, flux_timeseries = create_boundary_condition(condition_type, boundary,var_name,
                                  dxG, dyG, drF, hFacS, hFacW)

                if boundary=='north' and var_name=='VVEL':
                    flux_timeseries *=-1

                if not flux_started:
                    total_flux_timeseries = flux_timeseries
                    flux_started=True
                else:
                    total_flux_timeseries += flux_timeseries

                output_file = os.path.join(config_dir,'input','obcs_'+condition_type, 'BC_'+boundary+'_'+var_name+'_'+condition_type+'.bin')
                bc_grid.ravel('C').astype('>f4').tofile(output_file)

        cumulative_volume_timeseries = np.cumsum(3600*total_flux_timeseries)
        mean_etan_timeseries = cumulative_volume_timeseries/total_area

        ds = nc4.Dataset(os.path.join('..', 'input', 'expected_'+condition_type+'_results.nc'),'w')
        ds.createDimension('n',len(mean_etan_timeseries))
        var = ds.createVariable('volume','f4',('n',))
        var[:] = cumulative_volume_timeseries
        var = ds.createVariable('mean_etan', 'f4', ('n',))
        var[:] = mean_etan_timeseries
        ds.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_boundary_conditions(config_dir)
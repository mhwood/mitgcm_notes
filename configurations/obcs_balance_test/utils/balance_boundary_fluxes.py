
import os
import shutil
import simplegrid as sg
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
import argparse
import ast


def read_grid_information(config_dir):
    wetgrid_file = os.path.join(config_dir, 'input', 'grid.t001.nc')
    ds = nc4.Dataset(wetgrid_file)
    rA = np.array(ds.variables['rA'][:, :])
    dxG = np.array(ds.variables['dxG'][:, :])
    dyG = np.array(ds.variables['dyG'][:, :])
    hFacS = np.array(ds.variables['HFacS'][:, :, :])
    hFacW = np.array(ds.variables['HFacW'][:, :, :])
    delR = np.array(ds.variables['drF'][:])
    ds.close()
    return (rA, dxG, dyG, hFacS, hFacW, delR)

def read_velocities_normal_to_boundary(config_dir,experiment,Nr,nSx,nSy):

    unbalanced_obcs_folder = os.path.join(config_dir, 'input', 'obcs_'+experiment)

    uvel_west_file = os.path.join(unbalanced_obcs_folder,'BC_west_UVEL_'+experiment+'.bin')
    uvel_west = np.fromfile(uvel_west_file,'>f4')
    n_timesteps = int(np.size(uvel_west)/(Nr*nSy))
    uvel_west = np.reshape(uvel_west,(n_timesteps,Nr,nSy)).astype(float)

    # plt.subplot(1,3,1)
    # C = plt.imshow(uvel_west[17,:,:])
    # plt.colorbar(C)
    # plt.subplot(1, 3, 2)
    # C = plt.imshow(uvel_west[23, :, :])
    # plt.colorbar(C)
    # plt.subplot(1, 3, 3)
    # C = plt.imshow(uvel_west[23, :, :]!=0)
    # plt.colorbar(C)
    # plt.show()

    vvel_north_file = os.path.join(unbalanced_obcs_folder, 'BC_north_VVEL_'+experiment+'.bin')
    vvel_north = np.fromfile(vvel_north_file, '>f4')
    vvel_north = np.reshape(vvel_north, (n_timesteps, Nr, nSx)).astype(float)

    vvel_south_file = os.path.join(unbalanced_obcs_folder, 'BC_south_VVEL_'+experiment+'.bin')
    vvel_south = np.fromfile(vvel_south_file, '>f4')
    vvel_south = np.reshape(vvel_south, (n_timesteps, Nr, nSx)).astype(float)

    return(uvel_west,vvel_north,vvel_south)

def calculate_flux(boundary,timestep,field,width,delR,hFac,print_snapshot_stats):
    flux = np.zeros_like(field).astype(float)
    total_area = 0.0
    for i in range(np.shape(field)[1]):
        for j in range(np.shape(field)[0]):
            flux[j,i] = width[i]*delR[j]*hFac[j,i]*field[j,i]
            total_area += width[i]*delR[j]*hFac[j,i]

            if timestep==0 and print_snapshot_stats:
                if i==11 and boundary=='south' and j==82:
                    print('|- south -----------------------------|')
                    print('    drF(k) ',delR[j])
                    print('    hFacS[j,i] ',hFac[j,i])
                    print('    dxG[i] ',width[i])
                    print('    OBSv[j,i] ',field[j,i])

                if i==10 and boundary=='north' and j==68:
                    print('|- north -----------------------------|')
                    print('    drF(k) ',delR[j])
                    print('    hFacS[j,i] ',hFac[j,i])
                    print('    dxG[i] ',width[i])
                    print('    OBNv[j,i] ',field[j,i])

                if i==5 and boundary=='west' and j==85:
                    print('|- west -----------------------------|')
                    print('    drF(k) ',delR[j])
                    print('    hFacW[j,i] ',hFac[j,i])
                    print('    dyG[i] ',width[i])
                    print('    OBWu[j,i] ',field[j,i])

    total_flux = np.sum(flux)
    return(total_flux,total_area)

def calculate_flux_timeseries(uvel_west, vvel_north, vvel_south,
                              dxG, dyG, hFacS, hFacW, delR, print_snapshot_stats = False):

    boundaries = ['west','north','south']
    all_flux_timeseries = []
    all_areas = []
    all_masks = []
    for boundary in boundaries:
        if boundary=='north':
            hFac = hFacS[:,-2,:]
            width = dxG[-2,:]
            vel_grid = vvel_north
        if boundary=='south':
            hFac = hFacS[:,1,:]
            width = dxG[1, :]
            vel_grid = vvel_south
        if boundary=='west':
            hFac = hFacW[:,:,1]
            width = dyG[:, 1]
            vel_grid = uvel_west

        # note: the first and last grid cells are not calculated in the flux
        # one way to solve this issue is just to make hfac=0 in these areas
        hFac[:, 0] = 0
        hFac[:, -1] = 0

        # make a mask to apply the corrections later
        mask = np.copy(hFac)
        mask[mask>0]=1
        mask[mask<=0]=0
        all_masks.append(mask)

        n_timesteps = np.shape(vel_grid)[0]
        flux_timeseries = np.zeros((n_timesteps,)).astype(float)

        for timestep in range(n_timesteps):
            total_flux, total_area = calculate_flux(boundary, timestep, vel_grid[timestep,:,:], width, delR, hFac, print_snapshot_stats)
            flux_timeseries[timestep] = total_flux
        all_areas.append(total_area)

        if print_snapshot_stats:
            print('    Area ',total_area)
            print('    flow ',flux_timeseries[0])

        if boundary=='north':
            flux_timeseries *= -1

        all_flux_timeseries.append(flux_timeseries)

    return(all_flux_timeseries, all_areas, all_masks)

def plot_boundary_fluxes(output_file, all_flux_timeseries,rA):

    west_flux_timeseries = all_flux_timeseries[0]
    north_flux_timeseries = all_flux_timeseries[1]
    south_flux_timeseries = all_flux_timeseries[2]

    total_flux = west_flux_timeseries+north_flux_timeseries+south_flux_timeseries
    integrated_flux = np.cumsum(3600*total_flux)#/np.sum(rA)

    fig = plt.figure(figsize=(12, 14))
    plt.style.use("dark_background")

    plt.subplot(5, 1, 1)
    plt.plot(west_flux_timeseries/1e6)
    plt.title('West Boundary Flux Into Domain')
    plt.grid(linestyle='--', alpha=0.2)
    plt.gca().set_xticklabels([])
    plt.ylabel('Sv')

    plt.subplot(5, 1, 2)
    plt.plot(north_flux_timeseries/1e6)
    plt.title('North Boundary Flux Into Domain')
    plt.grid(linestyle='--', alpha=0.2)
    plt.gca().set_xticklabels([])
    plt.ylabel('Sv')

    plt.subplot(5, 1, 3)
    plt.plot(south_flux_timeseries/1e6)
    plt.title('South Boundary Flux Into Domain')
    plt.grid(linestyle='--', alpha=0.2)
    plt.gca().set_xticklabels([])
    plt.ylabel('Sv')

    plt.subplot(5, 1, 4)
    plt.plot(total_flux)
    plt.grid(linestyle='--', alpha=0.2)
    plt.title('Total Flux Into Domain')
    plt.gca().set_xticklabels([])
    plt.ylabel('Sv')

    plt.subplot(5, 1, 5)
    plt.plot(integrated_flux)
    plt.title('Integrated Flux Into Domain')
    plt.grid(linestyle='--', alpha=0.5)
    plt.ylabel('m (averaged over domain)')

    plt.savefig(output_file,bbox_inches='tight')
    plt.close(fig)


def balance_flux_on_boundaries(config_dir,experiment,running_average_radius,
                               uvel_west, vvel_north, vvel_south,
                               all_flux_timeseries, all_areas, all_masks):

    #####################################################
    # collect the info from each list

    west_flux_timeseries = all_flux_timeseries[0]
    north_flux_timeseries = all_flux_timeseries[1]
    south_flux_timeseries = all_flux_timeseries[2]

    west_area = all_areas[0]
    north_area = all_areas[1]
    south_area = all_areas[2]

    west_mask = all_masks[0]
    north_mask = all_masks[1]
    south_mask = all_masks[2]

    #####################################################
    # calculate the flux adjustment on each boundary
    # following the obcs package, a constant value is removed from all wet cells on the boundary
    # one can modify the balance fractions below, but the default it to remove it from all boundaries equally
    # note, for these, that the W and S boundaries have a -1 applied below

    balanceFacN = 1
    balanceFacS = 1
    balanceFacW = 1

    total_flux_timeseries = west_flux_timeseries + north_flux_timeseries + south_flux_timeseries
    total_area = west_area + north_area + south_area
    total_correction = total_flux_timeseries / total_area

    print('|- corrections -----------------------------|')
    print('    flowW ',west_flux_timeseries[0])
    print('    flowN ', north_flux_timeseries[0])
    print('    flowS ', south_flux_timeseries[0])
    print('    inFlow ',total_flux_timeseries[0])
    print('    areaOB ',total_area)
    print('    inFlow/areaOB ',total_correction[0])

    if running_average_radius>0:
        smooth_total_correction = np.zeros((len(total_correction),))
        for i in range(len(total_correction)):
            if i>=running_average_radius:
                smooth_total_correction[i] = np.mean(total_correction[i-running_average_radius:i])
        smooth_total_correction[:running_average_radius] = smooth_total_correction[running_average_radius]

        total_correction = smooth_total_correction

    # plt.plot(total_correction,'b-')
    # plt.plot(smooth_total_correction,'g-')
    # plt.show()

    #####################################################
    # apply the boundary flux adjustments to each boundary

    n_timesteps = np.size(total_flux_timeseries)

    print('    - Correcting the west flux')
    for i in range(n_timesteps):
        uvel_west[i,:,:] -= total_correction[i] * west_mask*balanceFacW
    output_file = os.path.join(config_dir, 'input', 'obcs_'+experiment[2:], 'BC_west_UVEL_'+experiment[2:]+'.bin')
    uvel_west.ravel('C').astype('>f4').tofile(output_file)

    print('    - Correcting the north flux')
    for i in range(n_timesteps):
        vvel_north[i, :, :] += total_correction[i] * north_mask * balanceFacN
    output_file = os.path.join(config_dir, 'input', 'obcs_' + experiment[2:], 'BC_north_VVEL_' + experiment[2:] + '.bin')
    vvel_north.ravel('C').astype('>f4').tofile(output_file)

    print('    - Correcting the south flux')
    for i in range(n_timesteps):
        vvel_south[i, :, :] -= total_correction[i] * south_mask * balanceFacS
    output_file = os.path.join(config_dir, 'input', 'obcs_' + experiment[2:],'BC_south_VVEL_' + experiment[2:] + '.bin')
    vvel_south.ravel('C').astype('>f4').tofile(output_file)

def balance_bc_fields(config_dir):

    Nr = 90
    nSx = 20
    nSy = 20

    experiment = 'unbalanced_constant'
    running_average_radius = 0

    # experiment = 'unbalanced_periodic'
    # running_average_radius = 24

    if 'obcs_'+experiment[2:] not in os.listdir(os.path.join(config_dir, 'input')):
        os.mkdir(os.path.join(config_dir, 'input', 'obcs_'+experiment[2:]))

    # these files are tangent to the boundary, so they dont need to be corrected
    if 'BC_north_UVEL_'+experiment[2:]+'.bin' not in os.listdir(os.path.join(config_dir, 'input', 'obcs_'+experiment[2:])):
        shutil.copyfile(os.path.join(config_dir, 'input', 'obcs_'+experiment,'BC_north_UVEL_'+experiment+'.bin'),
                   os.path.join(config_dir, 'input', 'obcs_'+experiment[2:],'BC_north_UVEL_'+experiment[2:]+'.bin'))
    if 'BC_south_UVEL_' + experiment[2:] + '.bin' not in os.listdir(
            os.path.join(config_dir, 'input', 'obcs_' + experiment[2:])):
        shutil.copyfile(os.path.join(config_dir, 'input', 'obcs_' + experiment, 'BC_south_UVEL_' + experiment + '.bin'),
                   os.path.join(config_dir, 'input', 'obcs_' + experiment[2:], 'BC_south_UVEL_' + experiment[2:] + '.bin'))
    if 'BC_west_VVEL_' + experiment[2:] + '.bin' not in os.listdir(
            os.path.join(config_dir, 'input', 'obcs_' + experiment[2:])):
        shutil.copyfile(os.path.join(config_dir, 'input', 'obcs_' + experiment, 'BC_west_VVEL_' + experiment + '.bin'),
                   os.path.join(config_dir, 'input', 'obcs_' + experiment[2:], 'BC_west_VVEL_' + experiment[2:] + '.bin'))

    print(' - Step 1: Reading in the grid information')
    rA, dxG, dyG, hFacS, hFacW, delR = read_grid_information(config_dir)

    print(' - Step 2: Reading in the unbalanced velocities normal to the boundaries')
    uvel_west, vvel_north, vvel_south = read_velocities_normal_to_boundary(config_dir,experiment,Nr,nSx,nSy)

    print(' - Step 3: Calculating the fluxes into the boundary')
    all_flux_timeseries, all_areas, all_masks  = \
        calculate_flux_timeseries(uvel_west, vvel_north, vvel_south,
                              dxG, dyG, hFacS, hFacW, delR, print_snapshot_stats = False)

    print(' - Step 4: Plotting the uncorrected fluxes')
    output_file = os.path.join(config_dir,'plots',experiment+'_fluxes.png')
    plot_boundary_fluxes(output_file, all_flux_timeseries,rA)

    print(' - Step 5: Adjusting the fluxes on each boundary')
    balance_flux_on_boundaries(config_dir, experiment, running_average_radius,
                               uvel_west, vvel_north, vvel_south,
                               all_flux_timeseries, all_areas, all_masks)

    print(' - Step 6: Reading in the now balanced velocities normal to the boundaries')
    uvel_west, vvel_north, vvel_south = read_velocities_normal_to_boundary(config_dir, experiment[2:], Nr, nSx, nSy)

    print(' - Step 7: Calculating the balanced fluxes into the boundary')
    all_flux_timeseries, all_areas, all_masks = \
        calculate_flux_timeseries(uvel_west, vvel_north, vvel_south,
                                  dxG, dyG, hFacS, hFacW, delR)

    print(' - Step 8: Plotting the balanced fluxes')
    output_file = os.path.join(config_dir, 'plots', experiment[2:] + '_fluxes.png')
    plot_boundary_fluxes(output_file, all_flux_timeseries, rA)



if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    balance_bc_fields(config_dir)


import os
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

def read_velocities_normal_to_boundary(config_dir,Nr,nSx,nSy):

    unbalanced_obcs_folder = os.path.join(config_dir, 'input', 'obcs_unbalanced')

    uvel_west_file = os.path.join(unbalanced_obcs_folder,'BC_west_UVEL_unbalanced.bin')
    uvel_west = np.fromfile(uvel_west_file,'>f4')
    n_timesteps = int(np.size(uvel_west)/(Nr*nSy))
    uvel_west = np.reshape(uvel_west,(n_timesteps,Nr,nSy))

    vvel_north_file = os.path.join(unbalanced_obcs_folder, 'BC_north_VVEL_unbalanced.bin')
    vvel_north = np.fromfile(vvel_north_file, '>f4')
    vvel_north = np.reshape(vvel_north, (n_timesteps, Nr, nSx))

    vvel_south_file = os.path.join(unbalanced_obcs_folder, 'BC_south_VVEL_unbalanced.bin')
    vvel_south = np.fromfile(vvel_south_file, '>f4')
    vvel_south = np.reshape(vvel_south, (n_timesteps, Nr, nSx))

    return(uvel_west,vvel_north,vvel_south)

def calculate_flux(field,width,delR,hFac):
    flux = np.zeros_like(field)
    total_area = 0
    for i in range(np.shape(field)[1]):
        for j in range(np.shape(field)[0]):
            flux[j,i] = width[i]*delR[j]*hFac[j,i]*field[j,i]
            total_area += width[i]*delR[j]*hFac[j,i]
    total_flux = np.sum(flux)
    return(total_flux,total_area)

def calculate_flux_timeseries(uvel_west, vvel_north, vvel_south,
                              dxG, dyG, hFacS, hFacW, delR):

    boundaries = ['west','north','south']
    all_flux_timeseries = []
    for boundary in boundaries:
        if boundary=='north':
            hFac = hFacS[:,-2,:]
            width = dxG[-2,:]
            vel_grid = uvel_west
        if boundary=='south':
            hFac = hFacS[:,1,:]
            width = dxG[1, :]
            vel_grid = vvel_north
        if boundary=='west':
            hFac = hFacW[:,:,1]
            width = dyG[:, 1]
            vel_grid = vvel_south

        # note: the first and last grid cells are not calculated in the flux
        # one way to solve this issue is just to make hfac=0 in these areas
        hFac[:, 0] = 0
        hFac[:, -1] = 0

        n_timesteps = np.shape(vel_grid)[0]
        flux_timeseries = np.zeros((n_timesteps,))

        for timestep in range(n_timesteps):
            total_flux, total_area = calculate_flux(vel_grid[timestep,:,:], width, delR, hFac)
            flux_timeseries[:] = total_flux

        if boundary=='north':
            flux_timeseries *= -1

        all_flux_timeseries.append(flux_timeseries)

    west_flux_timeseries = all_flux_timeseries[0]
    north_flux_timeseries = all_flux_timeseries[1]
    south_flux_timeseries = all_flux_timeseries[2]

    return(west_flux_timeseries, north_flux_timeseries, south_flux_timeseries)

def calculate_wet_boundary_area(boundary_name, bathy_grid, L1_wet_grid, DXC, DYC, delR):
    Z_bottom = np.cumsum(delR)
    Z_top = np.concatenate([np.array([0]), Z_bottom[:-1]])

    total_area = 0

    if boundary_name=='north':
        for i in range(np.shape(L1_wet_grid)[2]): # go along the boundary
            for d in range(len(Z_top)):
                if L1_wet_grid[d,-1,i]==1:
                    if bathy_grid[-1,i]<Z_bottom[d]:
                        frac=1
                    elif bathy_grid[-1,i]>Z_top[d] and bathy_grid[-1,i]<Z_bottom[d]:
                        frac = (bathy_grid[-1,i]-Z_top[d])/delR[d]
                    else:
                        frac=0
                    total_area+=frac*DXC[-1,i]*delR[d]

        area_grid = np.zeros((len(Z_top),np.shape(L1_wet_grid)[2]))
        for i in range(np.shape(L1_wet_grid)[2]): # go along the boundary
            for d in range(len(Z_top)):
                if L1_wet_grid[d,-1,i]==1:
                    if bathy_grid[-1,i]<Z_bottom[d]:
                        frac=1
                    elif bathy_grid[-1,i]>Z_top[d] and bathy_grid[-1,i]<Z_bottom[d]:
                        frac = (bathy_grid[-1,i]-Z_top[d])/delR[d]
                    else:
                        frac=0
                    area_grid[d,i] = frac*DXC[-1,i]*delR[d]


    if boundary_name=='south':
        for i in range(np.shape(L1_wet_grid)[2]): # go along the boundary
            for d in range(len(Z_top)):
                if L1_wet_grid[d,0,i]==1:
                    if bathy_grid[0,i]<Z_bottom[d]:
                        frac=1
                    elif bathy_grid[0,i]>Z_top[d] and bathy_grid[0,i]<Z_bottom[d]:
                        frac = (bathy_grid[0,i]-Z_top[d])/delR[d]
                    else:
                        frac=0
                    total_area+=frac*DXC[0,i]*delR[d]

        area_grid = np.zeros((len(Z_top), np.shape(L1_wet_grid)[2]))
        for i in range(np.shape(L1_wet_grid)[2]):  # go along the boundary
            for d in range(len(Z_top)):
                if L1_wet_grid[d, 0,i] == 1:
                    if bathy_grid[0,i] < Z_bottom[d]:
                        frac = 1
                    elif bathy_grid[0,i] > Z_top[d] and bathy_grid[0,i] < Z_bottom[d]:
                        frac = (bathy_grid[0,i] - Z_top[d]) / delR[d]
                    else:
                        frac = 0
                    area_grid[d, i] = frac * DXC[0,i] * delR[d]

    if boundary_name=='west':
        for i in range(np.shape(L1_wet_grid)[1]): # go along the boundary
            for d in range(len(Z_top)):
                if L1_wet_grid[d,i,0]==1:
                    if bathy_grid[i,0]<Z_bottom[d]:
                        frac=1
                    elif bathy_grid[i,0]>Z_top[d] and bathy_grid[i,0]<Z_bottom[d]:
                        frac = (bathy_grid[i,0]-Z_top[d])/delR[d]
                    else:
                        frac=0
                    total_area+=frac*DYC[i,0]*delR[d]

        area_grid = np.zeros((len(Z_top), np.shape(L1_wet_grid)[1]))
        for i in range(np.shape(L1_wet_grid)[1]):  # go along the boundary
            for d in range(len(Z_top)):
                if L1_wet_grid[d, i,0] == 1:
                    if bathy_grid[i,0] < Z_bottom[d]:
                        frac = 1
                    elif bathy_grid[i,0] > Z_top[d] and bathy_grid[i,0] < Z_bottom[d]:
                        frac = (bathy_grid[i,0] - Z_top[d]) / delR[d]
                    else:
                        frac = 0
                    area_grid[d, i] = frac * DYC[i,0] * delR[d]

    return(total_area, area_grid)

def plot_boundary_fluxes(output_file, west_flux_timeseries, north_flux_timeseries, south_flux_timeseries):

    total_flux = west_flux_timeseries+north_flux_timeseries+south_flux_timeseries
    integrated_flux = np.cumsum(total_flux)

    fig = plt.figure(figsize=(12, 14))
    plt.style.use("dark_background")

    plt.subplot(5, 1, 1)
    plt.plot(west_flux_timeseries,label='15 min interval',linewidth=0.5)
    plt.legend(ncol=2)
    plt.title('West Boundary Flux Into Domain')
    plt.grid(linestyle='--', alpha=0.5)
    plt.gca().set_xticklabels([])
    plt.ylabel('Sv')

    plt.subplot(5, 1, 2)
    plt.plot(north_flux_timeseries, label='15 min interval', linewidth=0.5)
    # plt.legend(ncol=2)
    plt.title('North Boundary Flux Into Domain')
    plt.grid(linestyle='--', alpha=0.5)
    plt.gca().set_xticklabels([])
    plt.ylabel('Sv')

    plt.subplot(5, 1, 3)
    plt.plot(south_flux_timeseries, label='15 min interval', linewidth=0.5)
    # plt.legend(ncol=2)
    plt.title('South Boundary Flux Into Domain')
    plt.grid(linestyle='--', alpha=0.5)
    plt.gca().set_xticklabels([])
    plt.ylabel('Sv')

    plt.subplot(5, 1, 4)
    plt.plot(total_flux, label='smoothed daily', linewidth=0.5)
    # plt.legend(ncol=2)
    plt.grid(linestyle='--', alpha=0.5)
    plt.gca().set_xticklabels([])
    plt.ylabel('Sv')

    plt.subplot(5, 1, 5)
    plt.plot(integrated_flux, label='raw', linewidth=2)
    plt.legend(ncol=2)
    plt.title('Integrated Flux Into Domain')
    plt.grid(linestyle='--', alpha=0.5)
    plt.ylabel('m (averaged over domain)')

    plt.savefig(output_file,bbox_inches='tight')
    plt.close(fig)

def balance_flux_on_boundaries(config_dir,n_rows_L1,n_cols_L1,delR,
                               total_smooth_flux,
                               north_fraction,south_fraction,west_fraction,
                               north_wet_area,south_wet_area,west_wet_area,
                               north_wet_grid,south_wet_grid,west_wet_grid):

    n_timesteps = np.size(total_smooth_flux)
    Nr = len(delR)

    print('    - Correcting the west flux')
    boundary_file = os.path.join(config_dir, 'L1_1080', 'input', 'obcs', 'L1_BC_west_UVEL_unbalanced.bin')
    boundary_grid = np.fromfile(boundary_file, '>f4')
    boundary_grid = np.reshape(boundary_grid, (n_timesteps, Nr, n_rows_L1))
    for i in range(n_timesteps):
        boundary_grid[i] -= west_fraction*(total_smooth_flux[i]/west_wet_area)*west_wet_grid
    output_file = os.path.join(config_dir, 'L1_1080', 'input', 'obcs', 'L1_BC_west_UVEL.bin')
    boundary_grid.ravel('C').astype('>f4').tofile(output_file)

    print('    - Correcting the north flux')
    boundary_file = os.path.join(config_dir, 'L1_1080', 'input', 'obcs', 'L1_BC_north_VVEL_unbalanced.bin')
    boundary_grid = np.fromfile(boundary_file, '>f4')
    boundary_grid = np.reshape(boundary_grid, (n_timesteps, Nr, n_cols_L1))
    for i in range(n_timesteps):
        boundary_grid[i] -= north_fraction*(total_smooth_flux[i]/north_wet_area)*north_wet_grid
    output_file = os.path.join(config_dir, 'L1_1080', 'input', 'obcs', 'L1_BC_north_VVEL.bin')
    boundary_grid.ravel('C').astype('>f4').tofile(output_file)

    print('    - Correcting the south flux')
    boundary_file = os.path.join(config_dir, 'L1_1080', 'input', 'obcs', 'L1_BC_south_VVEL_unbalanced.bin')
    boundary_grid = np.fromfile(boundary_file, '>f4')
    boundary_grid = np.reshape(boundary_grid, (n_timesteps, Nr, n_cols_L1))
    for i in range(n_timesteps):
        boundary_grid[i] -= south_fraction*(total_smooth_flux[i]/south_wet_area)*south_wet_grid
    output_file = os.path.join(config_dir, 'L1_1080', 'input', 'obcs', 'L1_BC_south_VVEL.bin')
    boundary_grid.ravel('C').astype('>f4').tofile(output_file)

def balance_bc_fields(config_dir):

    Nr = 90
    nSx = 20
    nSy = 20

    print(' - Step 0: Reading in the grid information')
    rA, dxG, dyG, hFacS, hFacW, delR = read_grid_information(config_dir)

    print(' - Step 1: Reading in the velocities normal to the boundaries')
    uvel_west, vvel_north, vvel_south = read_velocities_normal_to_boundary(config_dir,Nr,nSx,nSy)

    print(' - Step 2: Calculating the fluxes into the boundary')
    west_flux_timeseries, north_flux_timeseries, south_flux_timeseries  = \
        calculate_flux_timeseries(uvel_west, vvel_north, vvel_south,
                              dxG, dyG, hFacS, hFacW, delR)

    print(' - Step 3: Plotting the uncorrected fluxes')
    output_file = os.path.join(config_dir,'plots','unbalanced_fluxes.png')
    plot_boundary_fluxes(output_file, west_flux_timeseries, north_flux_timeseries, south_flux_timeseries)

    print(' - Step 3: Correcting the fluxes on each boundary ')

    print('     - Fractional areas of boundaries:')
    # print('       + North: ' + '{:.2f}'.format(north_fraction * 100) + ' %')
    # print('       + South: ' + '{:.2f}'.format(south_fraction * 100) + ' %')
    # print('       + West: ' + '{:.2f}'.format(west_fraction * 100) + ' %')

    print(' - Step 4: Reading in the velocities normal to the boundaries')

    print(' - Step 5: Plotting the corrected fluxes')



if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    balance_bc_fields(config_dir)

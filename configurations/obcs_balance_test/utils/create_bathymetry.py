
import os
import simplegrid as sg
import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import argparse

def map_the_bathymetry_to_subdomain(bathy_grid_faces,XC_faces,YC_faces,XC_subset,YC_subset,llc):

    bathy_45 = np.concatenate((bathy_grid_faces[4].T, bathy_grid_faces[5].T), axis=1)
    XC_45 = np.concatenate((XC_faces[4].T, XC_faces[5].T), axis=1)
    YC_45 = np.concatenate((YC_faces[4].T, YC_faces[5].T), axis=1)

    points = np.hstack([np.reshape(XC_45,(np.size(XC_45),1)),
                        np.reshape(YC_45,(np.size(YC_45),1))])
    values = np.reshape(bathy_45, (np.size(bathy_45), 1))

    grid = griddata(points,values,(XC_subset,YC_subset),method='nearest')

    return(grid)




########################################################################################################################

def create_bathymetry_file(ecco_dir):

    print('Creating the L1 bathymetry file')
    llc = 90

    n_rows = 20
    n_cols = 20

    # read in the LLC grid subset to faces
    grid_file_dir = os.path.join(ecco_dir, 'LLC' + str(llc) + '_Files', 'mitgrid_tiles')
    XC_faces = {}
    YC_faces = {}
    for i in range(1, 6):
        if i < 3:
            grid_dict = sg.gridio.read_mitgridfile(os.path.join(grid_file_dir, 'tile00' + str(i) + '.mitgrid'), llc,
                                                   3 * llc)
        if i == 3:
            grid_dict = sg.gridio.read_mitgridfile(os.path.join(grid_file_dir, 'tile00' + str(i) + '.mitgrid'), llc,
                                                   llc)
        if i > 3:
            grid_dict = sg.gridio.read_mitgridfile(os.path.join(grid_file_dir, 'tile00' + str(i) + '.mitgrid'), 3 * llc,
                                                   llc)
        XC_face = grid_dict['XC'].T
        YC_face = grid_dict['YC'].T
        XC_faces[i] = XC_face
        YC_faces[i] = YC_face

    # read in the LLC grid subset to faces
    bathy_file = os.path.join(ecco_dir,'LLC'+str(llc)+'_Files', 'input_init', 'bathy_llc'+str(llc))
    bathy_grid = np.fromfile(bathy_file,'>f4')
    bathy_grid_faces = {}
    points_counted = 0
    for i in range(1, 6):
        if i < 3:
            grid_subset = bathy_grid[points_counted:points_counted+3*llc*llc]
            bathy_grid_faces[i] = np.reshape(grid_subset,(3*llc,llc))
            points_counted += 3*llc*llc
        if i == 3 or i==6:
            grid_subset = bathy_grid[points_counted:points_counted + llc * llc]
            bathy_grid_faces[i] = np.reshape(grid_subset,(llc, llc))
            points_counted += llc * llc
        if i > 3 and i < 6:
            grid_subset = bathy_grid[points_counted:points_counted + 3 * llc * llc]
            bathy_grid_faces[i] = np.reshape(grid_subset,(llc, llc *3))
            points_counted += 3 * llc * llc

    # read in the grid subset to interpolate onto
    mitgrid_file = os.path.join('..','input','tile001.mitgrid')
    grid_dict = sg.gridio.read_mitgridfile(mitgrid_file,n_cols,n_rows)
    XC_subset = grid_dict['XC'].T
    YC_subset = grid_dict['YC'].T

    # # map the LLC files into the subdomain
    print('    Mapping the LLC'+str(llc)+' bathymetry onto the L1 domain')
    output_dir = os.path.join('..','input')
    grid = map_the_bathymetry_to_subdomain(bathy_grid_faces,XC_faces,YC_faces,XC_subset,YC_subset,llc)

    plt.imshow(grid,origin='lower',cmap='Blues_r')
    plt.show()

    output_file = os.path.join(output_dir, 'bathymetry.bin')
    grid.ravel('C').astype('>f4').tofile(output_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-e", "--ecco_directory", action="store",
                        help="Path to the ECCO directory where LLC files are stored.", dest="ecco_path",
                        type=str, required=True)

    args = parser.parse_args()
    ecco_path = args.ecco_path

    create_bathymetry_file(ecco_path)
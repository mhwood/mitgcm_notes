
import os
import simplegrid as sg
import numpy as np
import argparse


def create_boundary_condition(boundary, Nr, sNx, sNy):
    n_timestep = 2

    if boundary in ['north','south']:
        grid = np.zeros((n_timestep,Nr,sNx))
    else:
        grid = np.zeros((n_timestep, Nr, sNy))

    return(grid)


########################################################################################################################

def create_boundary_conditions(config_dir):

    print('Creating a set of temporary boundary conditions')
    Nr = 90
    sNx = 20
    sNy = 20

    condition_type = 'zero_flux'

    if 'obcs_'+condition_type not in os.listdir(os.path.join(config_dir,'input')):
        os.mkdir(os.path.join(config_dir,'input','obcs_'+condition_type))

    for boundary in ['west','north','south']:
        for var_name in ['UVEL','VVEL']:
            bc_grid = create_boundary_condition(boundary, Nr, sNx, sNy)

            output_file = os.path.join(config_dir,'input','obcs_'+condition_type, 'BC_'+boundary+'_'+var_name+'_'+condition_type+'.bin')
            bc_grid.ravel('C').astype('>f4').tofile(output_file)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the configuration is stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_boundary_conditions(config_dir)
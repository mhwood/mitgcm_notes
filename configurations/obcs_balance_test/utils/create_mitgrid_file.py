
import os
import simplegrid as sg
import numpy as np
import matplotlib.pyplot as plt
from ecco_v4_py import read_llc_to_tiles, read_llc_to_faces
import pyproj
import argparse


def read_in_L0_domain_components(L0_mitgrid_dir,llc):

    t4 = L0_mitgrid_dir + '/tile004.mitgrid'
    t5 = L0_mitgrid_dir + '/tile005.mitgrid'

    mg4 = sg.gridio.read_mitgridfile(t4, llc * 3, llc)
    mg5 = sg.gridio.read_mitgridfile(t5, llc * 3, llc)

    XC_45 = np.concatenate((mg4['XC'], mg5['XC']), axis=1)
    YC_45 = np.concatenate((mg4['YC'], mg5['YC']), axis=1)

    XG_45 = np.concatenate((mg4['XG'][:,:-1], mg5['XG']), axis=1)
    YG_45 = np.concatenate((mg4['YG'][:,:-1], mg5['YG']), axis=1)

    return(XC_45,YC_45,XG_45,YG_45)

def generate_L1_domain_bounds(XC_45,YC_45,XG_45,YG_45):

    # strategy here to is to expand out from SE corner,
    # approximately centering the domain on the SWOT calval point
    #calval_coords = (35.7, -125.4) # for reference

    SE1 = (32.5, -117)

    G = pyproj.Geod(ellps='sphere', a=1, b=1)
    ci_se1, cj_se1, cd_se1 = sg.util.nearest(SE1[1], SE1[0], XG_45, YG_45, G)

    # define NW corner tracer cell by an integer multiple
    # number of tiles from SE corner tracer cell,
    # approximately centered around calval point
    ci_nw1 = ci_se1 - 20 + 1
    cj_nw1 = cj_se1 - 20 + 1

    XG_L0_grid_L1_domain = XG_45[ci_nw1:ci_se1 + 2, cj_nw1:cj_se1 + 2]
    YG_L0_grid_L1_domain = YG_45[ci_nw1:ci_se1 + 2, cj_nw1:cj_se1 + 2]
    XC_L0_grid_L1_domain = XC_45[ci_nw1:ci_se1 + 1, cj_nw1:cj_se1 + 1]
    YC_L0_grid_L1_domain = YC_45[ci_nw1:ci_se1 + 1, cj_nw1:cj_se1 + 1]

    domain_coords = [XC_L0_grid_L1_domain,YC_L0_grid_L1_domain,XG_L0_grid_L1_domain,YG_L0_grid_L1_domain]

    return(domain_coords)

def generate_L2_domain_bounds_from_L1(L1_mitgrid):

    # strategy here to is to expand out from SE corner,
    # approximately centering the domain on the SWOT calval point
    #calval_coords = (35.7, -125.4) # for reference

    SE2 = (30.4, -115)
    # SE2 = (27.4, -113)

    XC_L1 = L1_mitgrid['XC'].T
    YC_L1 = L1_mitgrid['YC'].T
    XG_L1 = L1_mitgrid['XG'].T
    YG_L1 = L1_mitgrid['YG'].T

    G = pyproj.Geod(ellps='sphere', a=1, b=1)
    ci_se2, cj_se2, cd_se2 = sg.util.nearest(SE2[1], SE2[0], XG_L1, YG_L1, G)

    # define NW corner tracer cell by an integer multiple
    # number of tiles from SE corner tracer cell,
    # approximately centered around calval point
    ci_nw2 = ci_se2 + 180 - 1
    cj_nw2 = cj_se2 - 240 + 1

    XG_L1_grid_L2_domain = XG_L1[ci_se2:ci_nw2 + 2, cj_nw2:cj_se2 + 2]
    YG_L1_grid_L2_domain = YG_L1[ci_se2:ci_nw2 + 2, cj_nw2:cj_se2 + 2]
    XC_L1_grid_L2_domain = XC_L1[ci_se2:ci_nw2 + 1, cj_nw2:cj_se2 + 1]
    YC_L1_grid_L2_domain = YC_L1[ci_se2:ci_nw2 + 1, cj_nw2:cj_se2 + 1]

    # plt.plot(XC_L1_grid_L2_domain,YC_L1_grid_L2_domain,'g.',markersize=6)
    # plt.plot(XC_L1,YC_L1,'k.',markersize=3)
    # plt.show()

    # flip everything so it turns out the right way up
    XG_L1_grid_L2_domain = np.flipud(XG_L1_grid_L2_domain)
    YG_L1_grid_L2_domain = np.flipud(YG_L1_grid_L2_domain)
    XC_L1_grid_L2_domain = np.flipud(XC_L1_grid_L2_domain)
    YC_L1_grid_L2_domain = np.flipud(YC_L1_grid_L2_domain)

    domain_coords = [XC_L1_grid_L2_domain, YC_L1_grid_L2_domain, XG_L1_grid_L2_domain, YG_L1_grid_L2_domain]

    # now, extend the domain coords - these will be used to create diagnostics_vec_masks

    ci_nw2 += 1
    ci_se2 -= 1
    cj_nw2 -= 1
    cj_se2 += 1

    XG_L1_grid_L2_domain_extended = XG_L1[ci_se2:ci_nw2 + 2, cj_nw2:cj_se2 + 2]
    YG_L1_grid_L2_domain_extended = YG_L1[ci_se2:ci_nw2 + 2, cj_nw2:cj_se2 + 2]
    XC_L1_grid_L2_domain_extended = XC_L1[ci_se2:ci_nw2 + 1, cj_nw2:cj_se2 + 1]
    YC_L1_grid_L2_domain_extended = YC_L1[ci_se2:ci_nw2 + 1, cj_nw2:cj_se2 + 1]

    # flip everything so it turns out the right way up
    XG_L1_grid_L2_domain_extended = np.flipud(XG_L1_grid_L2_domain_extended)
    YG_L1_grid_L2_domain_extended = np.flipud(YG_L1_grid_L2_domain_extended)
    XC_L1_grid_L2_domain_extended = np.flipud(XC_L1_grid_L2_domain_extended)
    YC_L1_grid_L2_domain_extended = np.flipud(YC_L1_grid_L2_domain_extended)

    domain_coords_extended = [XC_L1_grid_L2_domain_extended, YC_L1_grid_L2_domain_extended,
                              XG_L1_grid_L2_domain_extended, YG_L1_grid_L2_domain_extended]
    # domain_coords = domain_coords_extended = []
    return (domain_coords, domain_coords_extended)

def generate_L3_domain_bounds_from_L2(L2_mitgrid):

    # strategy here to is to expand out from SE corner,
    # approximately centering the domain on the SWOT calval point
    #calval_coords = (35.7, -125.4) # for reference

    SE2 = (32.5, -117)
    # SE2 = (27.4, -113)

    XC_L2 = L2_mitgrid['XC'].T
    YC_L2 = L2_mitgrid['YC'].T
    XG_L2 = L2_mitgrid['XG'].T
    YG_L2 = L2_mitgrid['YG'].T

    G = pyproj.Geod(ellps='sphere', a=1, b=1)
    ci_se2, cj_se2, cd_se2 = sg.util.nearest(SE2[1], SE2[0], XG_L2, YG_L2, G)

    x_shift = 0
    y_shift = 0
    cj_se2 += x_shift

    # define NW corner tracer cell by an integer multiple
    # number of tiles from SE corner tracer cell,
    # approximately centered around calval point
    ci_nw2 = ci_se2 + 240 - 1 + y_shift
    cj_nw2 = cj_se2 - 300- x_shift + 1

    XG_L2_grid_L3_domain = XG_L2[ci_se2:ci_nw2 + 2, cj_nw2:cj_se2 + 2]
    YG_L2_grid_L3_domain = YG_L2[ci_se2:ci_nw2 + 2, cj_nw2:cj_se2 + 2]
    XC_L2_grid_L3_domain = XC_L2[ci_se2:ci_nw2 + 1, cj_nw2:cj_se2 + 1]
    YC_L2_grid_L3_domain = YC_L2[ci_se2:ci_nw2 + 1, cj_nw2:cj_se2 + 1]

    # plt.plot(XC_L2_grid_L3_domain,YC_L2_grid_L3_domain,'g.',markersize=6)
    # plt.plot(XC_L2,YC_L2,'k.',markersize=3)
    # plt.show()

    # flip everything so it turns out the right way up
    XG_L2_grid_L3_domain = np.flipud(XG_L2_grid_L3_domain)
    YG_L2_grid_L3_domain = np.flipud(YG_L2_grid_L3_domain)
    XC_L2_grid_L3_domain = np.flipud(XC_L2_grid_L3_domain)
    YC_L2_grid_L3_domain = np.flipud(YC_L2_grid_L3_domain)

    domain_coords = [XC_L2_grid_L3_domain, YC_L2_grid_L3_domain, XG_L2_grid_L3_domain, YG_L2_grid_L3_domain]

    # now, extend the domain coords - these will be used to create diagnostics_vec_masks

    ci_nw2 += 1
    ci_se2 -= 1
    cj_nw2 -= 1
    cj_se2 += 1

    XG_L2_grid_L3_domain_extended = XG_L2[ci_se2:ci_nw2 + 2, cj_nw2:cj_se2 + 2]
    YG_L2_grid_L3_domain_extended = YG_L2[ci_se2:ci_nw2 + 2, cj_nw2:cj_se2 + 2]
    XC_L2_grid_L3_domain_extended = XC_L2[ci_se2:ci_nw2 + 1, cj_nw2:cj_se2 + 1]
    YC_L2_grid_L3_domain_extended = YC_L2[ci_se2:ci_nw2 + 1, cj_nw2:cj_se2 + 1]

    # flip everything so it turns out the right way up
    XG_L2_grid_L3_domain_extended = np.flipud(XG_L2_grid_L3_domain_extended)
    YG_L2_grid_L3_domain_extended = np.flipud(YG_L2_grid_L3_domain_extended)
    XC_L2_grid_L3_domain_extended = np.flipud(XC_L2_grid_L3_domain_extended)
    YC_L2_grid_L3_domain_extended = np.flipud(YC_L2_grid_L3_domain_extended)

    domain_coords_extended = [XC_L2_grid_L3_domain_extended, YC_L2_grid_L3_domain_extended,
                              XG_L2_grid_L3_domain_extended, YG_L2_grid_L3_domain_extended]
    # domain_coords = domain_coords_extended = []
    return (domain_coords, domain_coords_extended)


def create_mitgrid_matrices(domain_coords):
    mitgrid_matrices = dict()
    mitgrid_matrices['XG'] = domain_coords[2]
    mitgrid_matrices['YG'] = domain_coords[3]
    mitgrid_matrices['XC'] = domain_coords[0]
    mitgrid_matrices['YC'] = domain_coords[1]
    return(mitgrid_matrices)

def create_new_mitgrid(output_file,mitgrid_matrices_L0_grid_L1_domain,XG_L0_grid_L1_domain,YG_L0_grid_L1_domain,factor):
    mg_new_L0_grid_L1_domain, n_rows, n_cols = \
        sg.regrid.regrid(mitgrid_matrices=mitgrid_matrices_L0_grid_L1_domain, \
                         lon_subscale=factor, lat_subscale=factor, \
                         lon1=XG_L0_grid_L1_domain[0, 0], lat1=YG_L0_grid_L1_domain[0, 0],
                         lon2=XG_L0_grid_L1_domain[-1, -1], lat2=YG_L0_grid_L1_domain[-1, -1],
                         verbose=False, outfile = output_file)

    sg.gridio.write_mitgridfile(output_file, mg_new_L0_grid_L1_domain, n_rows, n_cols)
    return(mg_new_L0_grid_L1_domain, n_rows, n_cols)

def create_mitgrid_files(ecco_path,print_status,L0_llc):

    shape_dict = {}

    # make a list of pertinent paths
    L0_llc_mitgrid_dir = os.path.join(ecco_path,'LLC'+str(L0_llc)+'_Files','mitgrid_tiles')

    # step 1: read in bathy faces overlapping region of interest
    if print_status:
        print('Reading faces 4 and 5 from the LLC'+str(L0_llc)+' (L0) configuration')
    XC_45,YC_45,XG_45,YG_45 = read_in_L0_domain_components(L0_llc_mitgrid_dir,L0_llc)

    # step 2: expand the L1 domain
    L1_domain_coords = generate_L1_domain_bounds(XC_45,YC_45,XG_45,YG_45)

    mitgrid_matrices_L0_grid_L1_domain = create_mitgrid_matrices(L1_domain_coords)

    output_dir = '..'
    if print_status:
        print('    Generating mitgrid for the L1_'+str(L0_llc)+' domain')
    output_file = os.path.join(output_dir,'input','tile001.mitgrid')
    L0_llc_mitgrid, rows, cols = create_new_mitgrid(output_file, mitgrid_matrices_L0_grid_L1_domain,
                       L1_domain_coords[2],L1_domain_coords[3],factor=1)
    if print_status:
        print('        The L1_'+str(L0_llc)+' domain has '+str(cols)+' rows and '+str(rows)+' cols')



if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-e", "--ecco_directory", action="store",
                        help="Path to the ECCO directory where LLC files are stored.", dest="ecco_path",
                        type=str, required=False, default = '../ECCO')

    parser.add_argument("-l", "--llc", action="store",
                        help="The LLC number of the source domain.", dest="llc",
                        type=int, required=True)

    parser.add_argument("-p", "--print_status", action="store",
                        help="Print status of routine (1 for True, 0 for False).", dest="print_status",
                        type=int, required=False, default=1)

    args = parser.parse_args()
    ecco_path = args.ecco_path
    L0_llc = args.llc
    print_status = args.print_status

    if print_status>0:
        print_status=True
    else:
        print_status=False

    create_mitgrid_files(ecco_path,print_status,L0_llc)

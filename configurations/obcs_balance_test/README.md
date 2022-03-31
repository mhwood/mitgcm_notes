# obcs_balance_test

## Purpose
This configuration was created to test offline Python code to balance fluxes into a regional domain. This was motivated by a need to remove a mean volume flux into a domain, rather than balancing at each time step. This functionality is not currently available in the obcs package.

## Description
This model is configured in a 20 x 20 model grid with 90 vertical levels. The bathymetry is derived from ECCO's LLC90 grid and there are no external forcings applied to the surface. The boundary conditions are varied in several different experiements (see below).

## Building the model and retrieving the grid
There are two parts to set up this model and prepare for the expeiment. First, some model files will be created and the model will be built. Next, the model will be run for 1 timestep to retrieve the model-generate grid parametets (particularly the hFac grids). The steps to run and assess the experiements is described below

### Step 1: Set up MITgcm
Add this `obcs_balance_test` directory to a (preferrably fresh) clone of the MITgcm, as `MITgcm/configurations/obcs_balance_test`

The following scripts rely on access to the `config_dir`. If you `cd` to the `obcs_balance_test/utils` directory, then the `configdir = ../`

### Step 2: Create the model grid 

This component is carried out with the simplegrid package:
```
python create_mitgrid_file.py -d config_dir -e ecco_dir
```
Note: this script was written with reference to a directory storing ECCO grid files. In particular, it uses files `til004.mitgrid` and `til005.mitgrid` located in `${ecco_dir}/LLC90_Files/mitgrid_tiles`

### Step 3: Create boundary conditions
Create some boundary conditions which are just 0's (remember we are mostly using this setup to dump out the grid)
```
python create_zero_boundary_conditions.py -d config_dir
```

### Step 4: Copy the grid to the input directory
```
cp ../run_zero_flux/mnc_0001/grid.t001.nc ../input
```

## Run the model experiments
Now that the model is built and the grid is obtained, boundary conditions can be created, and the balancing codes can be tested.


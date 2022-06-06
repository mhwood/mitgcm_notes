# precip_slr_corr

Cliff notes to clean up and test
- clone MITgcm
- put this dir in MITgcm/configurations/
- move slr package to pkg dir
- edit PARAMS.h to have the keywords
- copy tutorial data to the data dir
- cd to utils
- build using bash script

Run the spin-up model first
- run with the run script
    - link the input files to the run dir
    - link the data files to the run dir
- plot status with a python script if desired

Run the balanced model
- run with the run script
    - link the input files to the run dir
    - link the data files to the run dir
    - link the pickup file into the run dir
- plot status with a python script if desired


# precip_slr_corr

## Purpose
This configuration is designed to test a new package (```slr_corr```) which can scale precipitation fields within the uncertainties of the precipitation field in order to match observed sea level. 

## Description
This model is constructed by modifying the global_with_exf tutorial. The `ncep_qnet` and `ncep_emp` files are removed and replaced with evaporation and precipitation files constructed with ECCO data.

### Testing the package
To test this package, first the model is spun up for 20 years. Then, a number of configurations are available to test different balance intervals relative to a given sea level curve.

## Step 1: Set up MITgcm, build the model, run the spin-up
Add this `precip_slr_corr` directory to a (preferrably fresh) clone of the MITgcm, as `MITgcm/configurations/precip_slr_corr`

The following scripts rely on relative paths from inside the utils director. To build the model, use the following script:
```
cd MITgcm/configurations/precip_slr_corr/utils
bash build.sh
```
Then, collect/create the pertinent data files for the spin-up, and run the spin-up model using the following commands:
```
python3 create_steady_state_files.py
bash run_spinup.sh
```

## Step 2: Run the SLR scenario
Next, we create a precipitation file which will induce a 1 mm/yr sea level rise signal (with a seasonal cycle).



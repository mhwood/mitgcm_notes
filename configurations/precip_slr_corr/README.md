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
This model is configured is constructed by modifying the global_with_exf tutorial. The ncep_qnet and ncep_emp files are removed and replaced with evaporation and precipitation files constructed with ECCO data.

### Testing the package
To test this package, first the model is spun up for 20 years. Then, a number of configurations are available to test different balance intervals relative to a given sea level curve.

## Step 1: Set up MITgcm
Add this `precip_slr_corr` directory to a (preferrably fresh) clone of the MITgcm, as `MITgcm/configurations/precip_slr_corr`

The following scripts rely on relative paths from inside the utils director. To build the model, use the following script:
```
bash build.sh
```

# precip_balance_corr

This configuration is designed to test new lackage (```slr```) which can scale precipitation fields within the uncertainties of the precipitation field. 

Note to self: fix ECCO file so there is no negative precip, then remove useExfCheckRange=.FALSE. from data.exf

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

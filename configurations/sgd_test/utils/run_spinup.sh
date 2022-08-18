cd ../
mkdir run_spinup
rm -r run_spinup/mnc_0001
rm run_spinup/*
cd run_spinup
ln -s ../input/* .
ln -s ../namelist_spinup/* .
ln -s ../build/mitgcmuv .
mpirun -np 1 ./mitgcmuv

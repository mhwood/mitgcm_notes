cd ../
mkdir run_spinup
rm -r run_spinup/mnc_0001
rm -r run_spinup/mnc_0002
rm run_spinup/*
cd run_spinup
ln -s ../input_spinup/* .
ln -s ../data/* .
ln -s ../build/mitgcmuv .
mpirun -np 2 ./mitgcmuv

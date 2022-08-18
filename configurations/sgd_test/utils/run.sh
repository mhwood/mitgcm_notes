cd ../
mkdir run
rm -r run/mnc_0001
rm run/*
cd run
ln -s ../input/* .
ln -s ../namelist/* .
ln -s ../build/mitgcmuv .
mpirun -np 1 ./mitgcmuv > output.txt

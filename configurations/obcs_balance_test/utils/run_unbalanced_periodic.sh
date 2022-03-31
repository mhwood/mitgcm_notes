cd ../
mkdir run_unbalanced_periodic
rm -r run_unbalanced_periodic/mnc_0001
rm run_unbalanced_periodic/*
cd run_unbalanced_periodic
ln -s ../input/* .
mv data.obcs_unbalanced_periodic data.obcs
mv obcs_unbalanced_periodic obcs
ln -s ../build/mitgcmuv .
./mitgcmuv > output.txt
cd ../
mkdir run_balanced_periodic_manual
rm -r run_balanced_periodic_manual/mnc_0001
rm run_balanced_periodic_manual/*
cd run_balanced_periodic_manual
ln -s ../input/* .
mv data.obcs_balanced_periodic_manual data.obcs
mv obcs_balanced_periodic obcs
ln -s ../build/mitgcmuv .
./mitgcmuv > output.txt
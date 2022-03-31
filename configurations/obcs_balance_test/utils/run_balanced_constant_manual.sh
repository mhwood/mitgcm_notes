cd ../
mkdir run_balanced_constant_manual
rm -r run_balanced_constant_manual/mnc_0001
rm run_balanced_constant_manual/*
cd run_balanced_constant_manual
ln -s ../input/* .
mv data.obcs_balanced_constant_manual data.obcs
mv obcs_balanced_constant obcs
ln -s ../build/mitgcmuv .
./mitgcmuv > output.txt
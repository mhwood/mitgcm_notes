cd ../
mkdir run_balanced_constant_by_mitgcm
rm -r run_balanced_constant_by_mitgcm/mnc_0001
rm run_balanced_constant_by_mitgcm/*
cd run_balanced_constant_by_mitgcm
ln -s ../input/* .
mv data.obcs_balanced_constant_by_mitgcm data.obcs
mv obcs_unbalanced_constant obcs
ln -s ../build/mitgcmuv .
./mitgcmuv > output.txt
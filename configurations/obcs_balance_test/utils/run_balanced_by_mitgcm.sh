cd ../
rm -r run_balanced_by_mitgcm/mnc_0001
rm run_balanced_by_mitgcm/*
cd run_balanced_by_mitgcm
ln -s ../input/* .
mv data.obcs_balanced_by_mitgcm data.obcs
mv obcs_unbalanced obcs
ln -s ../build/mitgcmuv .
./mitgcmuv > output.txt
cd ../
rm -r run_unbalanced/mnc_0001
rm run_unbalanced/*
cd run_unbalanced
ln -s ../input/* .
mv data.obcs_unbalanced data.obcs
mv obcs_unbalanced obcs
ln -s ../build/mitgcmuv .
./mitgcmuv > output.txt
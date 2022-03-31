cd ../
mkdir run_unbalanced_constant
rm -r run_unbalanced_constant/mnc_0001
rm run_unbalanced_constant/*
cd run_unbalanced_constant
ln -s ../input/* .
mv data.obcs_unbalanced_constant data.obcs
mv obcs_unbalanced_constant obcs
ln -s ../build/mitgcmuv .
./mitgcmuv > output.txt
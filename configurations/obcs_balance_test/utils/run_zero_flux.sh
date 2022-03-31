cd ../
mkdir run_zero_flux
rm -r run_zero_flux/mnc_0001
rm run_zero_flux/*
cd run_zero_flux
ln -s ../input/* .
cp ../input/data_zero_flux data
mv data.obcs_zero_flux data.obcs
mv obcs_zero_flux obcs
ln -s ../build/mitgcmuv .
./mitgcmuv > output.txt
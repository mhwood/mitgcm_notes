

import os
import numpy as np
import netCDF4 as nc4


config_dir = '../'

temp_profile = []

n_rows = 31
n_cols = 90

temp_profile = np.array([ 5.22489767e+00, 1.24813611e+00, 1.35338919e-03, -7.22065829e-01,
                -1.28098373e+00, -1.54885511e+00, -1.58127672e+00, -1.54374410e+00,
                -1.48511965e+00, -1.36213107e+00, -1.31243087e+00, -1.31935274e+00,
                -1.04217728e+00, -8.90446470e-01, -8.53148577e-01, -7.56841498e-01,
                -7.24989401e-01, -6.69863164e-01, -5.75914957e-01, -4.72008116e-01,
                -2.77194279e-01, -7.35417611e-02, 4.23780354e-02, 1.76672296e-01,
                3.14941406e-01, 3.95308091e-01, 4.82657981e-01, 5.81283511e-01,
                6.75010176e-01, 7.16658061e-01, 7.66042318e-01, 8.25380487e-01,
                8.74738518e-01, 8.94170028e-01, 9.23568915e-01, 9.42996009e-01,
                9.41494416e-01, 9.32000002e-01, 9.27061328e-01, 9.11039559e-01,
                9.04992022e-01, 8.94511536e-01, 8.79610605e-01, 8.69135492e-01,
                8.64200022e-01, 8.58148520e-01, 8.47668239e-01, 8.47146440e-01,
                8.45625434e-01, 8.36143922e-01])

salt_profile = np.array([25.12191108, 29.81428471, 31.23700725, 32.01968108, 32.51934312, 32.86473642,
                 33.10287285, 33.2662957,  33.35362318, 33.46010635, 33.57328894, 33.6873636,
                 33.81224367, 33.92033494, 34.02783166, 34.14469759, 34.19523662, 34.26869987,
                 34.30926476, 34.35877634, 34.4305267,  34.49247861, 34.52588726, 34.57256092,
                 34.59973035, 34.64439124, 34.65740346, 34.69040878, 34.73338395, 34.73068791,
                 34.74809313, 34.77000406, 34.7944238,  34.81634208, 34.84146189, 34.85186144,
                 34.85222826, 34.87700785, 34.87604708, 34.88751185, 34.88780715, 34.8930785,
                 34.88899005, 34.88853945, 34.88188353, 34.88218082, 34.88748682, 34.89019355,
                 34.89488542, 34.89246761])


Nr = 50

temp_grid = np.zeros((Nr, n_rows, n_cols))
salt_grid = np.zeros((Nr, n_rows, n_cols))

for i in range(n_rows):
    for j in range(n_cols):
        temp_grid[:,i,j] = temp_profile
        salt_grid[:,i,j] = salt_profile

output_file = os.path.join(config_dir,'input','polar_theta_IC.bin')
temp_grid.ravel(order='C').astype('>f4').tofile(output_file)

output_file = os.path.join(config_dir,'input','polar_salinity_IC.bin')
salt_grid.ravel(order='C').astype('>f4').tofile(output_file)

C------------------------------------------------------------------------------|
C                           SLR_CORR_SIZE.h
C------------------------------------------------------------------------------|
C This determines the maximum size of the vector which records EtaN errors
C The first N values are averaged when SLRC_balancePeriod > 0 
C The length needs to be a minimum of ceil(SLRC_balancePeriod/timestep)

      INTEGER, PARAMETER :: SLRC_n_rec = 6720


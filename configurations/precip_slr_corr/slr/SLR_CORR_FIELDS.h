

CBOP
C     !ROUTINE: SLRC_CORR_FIELDS.h
C     !INTERFACE:
C     #include "SLR_CORR_FIELDS.h"

C     !DESCRIPTION:
C     *==========================================================*
C     | SLR_CORR_FIELDS.h
C     | o Header file containing EtaN errors
C     *==========================================================*
CEOP

#ifdef ALLOW_SLR_CORR
      COMMON /SLR_CORR_FIELDS/
     & SLRC_errors
      _RL SLRC_errors(SLRC_n_rec)
#endif /* ALLOW_SLR_CORR */

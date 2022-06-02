#ifdef ALLOW_SLR_CORR

CBOP
C     !ROUTINE: SLR_CORR_PARAMS.h
C     !INTERFACE:
C     #include "SLR_CORR_PARAMS.h"

C     !DESCRIPTION:
C     *==========================================================*
C     | SLR_CORR_PARAMS.h
C     | o Header file containing SLR_CORR parameters
C     *==========================================================*
C     | o Note: does not (and should not) contain any conditional
C     |   statement that depends on SLR_CORR options ; therefore
C     |   can be safely included without SLR_CORR_OPTIONS.h
C     *==========================================================*
CEOP

C--   COMMON /SLRC_PARM_L/ SLR_CORR logical-type parameter

      COMMON /SLRC_PARM_L/
     & useSLR_CORRbalance
      LOGICAL useSLR_CORRbalance

C--   COMMON /SLRC_PARM_R/ SLR_CORR real-type parameter
      COMMON /SLRC_PARM_R/
     &     SLR_CORR_balancePeriod
      _RL SLR_CORR_balancePeriod

C--   COMMON /SLRC_FILES/ SLR_CORR character-type parameter
      COMMON /SLRC_FILES/
     &      SLR_target_File
      CHARACTER*(MAX_LEN_FNAM)
     &      SLR_target_File

#endif /* ALLOW_SLR_CORR */

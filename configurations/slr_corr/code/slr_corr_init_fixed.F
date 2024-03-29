C----&------------------------------------------------------------------xxxxxxx|
C                           slr_corr_init_fixed
C----&------------------------------------------------------------------xxxxxxx|

#include "SLR_CORR_OPTIONS.h"

C----&------------------------------------------------------------------xxxxxxx|
C----&------------------------------------------------------------------xxxxxxx|
CBOP
C !ROUTINE: SLR_CORR_INIT_FIXED

C !INTERFACE:
      SUBROUTINE SLR_CORR_INIT_FIXED( myThid )

C     !DESCRIPTION:
C     Initialize SLR_CORR parameters and constants.

C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "SLR_CORR_PARAM.h"
#include "SLR_CORR_FIELDS.h"

C     !INPUT PARAMETERS:
      INTEGER myThid

C     !LOCAL VARIABLES:
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      INTEGER iUnit
      INTEGER, PARAMETER :: debug = 1

CEOP
C----&------------------------------------------------------------------xxxxxxx|
C----&------------------------------------------------------------------xxxxxxx|
CBOC

#ifdef ALLOW_SLR_CORR

C----&------------------------------------------------------------------xxxxxxx|
C     Define the slr_corr variables
C----&------------------------------------------------------------------xxxxxxx|

      INTEGER slrc_date_array(4)
      INTEGER slrc_difftime(4)
      INTEGER slrc_startdate(4)
      
C !MORE LOCAL VARIABLES:
      INTEGER iRec
      _RS     dummyRS(1)
      CHARACTER*(2) fType
      INTEGER ioUnit, bi, bj
      CHARACTER*(MAX_LEN_FNAM) fullName
    
      if (debug .eq. 1) then
        WRITE(msgBuf,'(A)') "========================================="
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
        WRITE(msgBuf,'(A)') "  BEGIN SLR_CORR_INIT_FIXED DEBUG "
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
        WRITE(msgBuf,'(A)') "========================================="
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
        endif

C----&------------------------------------------------------------------xxxxxxx|
C   Calulate the start time of the observations based on the input dates
C----&------------------------------------------------------------------xxxxxxx|

      _BEGIN_MASTER(myThid)

      CALL CAL_FULLDATE( slrc_obs_startdate_1, slrc_obs_startdate_2,
     &                     slrc_date_array, myThid )
      CALL cal_getdate( nIter0, startTime, slrc_startdate, myThid )
      CALL cal_TimePassed( slrc_startdate, slrc_date_array, 
     &                     slrc_difftime,
     I                     myThid )
      CALL cal_ToSeconds ( slrc_difftime, slrc_obs_start_time, myThid )
      slrc_obs_start_time = startTime  + slrc_obs_start_time

C      IF ( ( fld_startdate1.NE.0 .OR. fld_startdate2.NE.0 )
C     &      .AND. fld_period.GT.0. ) THEN
C-    Report inconsistent setting:
C         IF ( fld_startdate1.NE.0 ) THEN
C          WRITE(msgBuf,'(8A)') 'S/R EXF_GETFFIELD_START: ',
C     &   'start-date for ', pkg_name, '-field "', fld_name,
C     &   '" = ', fld_name, 'startdate1'
C          CALL PRINT_ERROR( msgBuf, myThid )
C         ENDIF
C         IF ( fld_startdate2.NE.0 ) THEN
C          WRITE(msgBuf,'(8A)') 'S/R EXF_GETFFIELD_START: ',
C     &   'start-date for ', pkg_name, '-field "', fld_name,
C     &   '" = ', fld_name, 'startdate2'
C          CALL PRINT_ERROR( msgBuf, myThid )
C         ENDIF
C         WRITE(msgBuf,'(5A)') '    ',
C     &   'is not allowed (in data.', pkg_name, ')',
C     &   ' when pkg/cal is not used (useCAL=F)'
C         CALL PRINT_ERROR( msgBuf, myThid )
C         errCount = errCount + 1
C        ENDIF

       PRINT *, 'slrc_obs_start_time',slrc_obs_start_time

C----&------------------------------------------------------------------xxxxxxx|
C   Read in the target ETAN timeseries (mean slr observations)
C----&------------------------------------------------------------------xxxxxxx|

      ioUnit = 0
      fType ='RL'
      bi = 0
      bj = 0
      iRec = 1

#ifdef ALLOW_MDSIO
      CALL MDS_READVEC_LOC(
     &                     slrc_obs_filename, slrc_filePrec, ioUnit,
     &                     fType, slrc_n_obs, 
     &                     slrc_obs_timeseries, dummyRS,
     I                     bi, bj, iRec, myThid )
#else
      STOP 'ABNORMAL END: S/R READ_GLVEC_RL needs MDSIO pkg'
#endif

      PRINT *, 'slrc_obs_timeseries',slrc_obs_timeseries

C----&------------------------------------------------------------------xxxxxxx|
C     Wrap up the code
C----&------------------------------------------------------------------xxxxxxx|

      _END_MASTER(myThid)

C--   Everyone must wait for the parameters to be loaded
      _BARRIER

#endif /* ALLOW_SLR_CORR */

      if (debug .eq. 1) then
        WRITE(msgBuf,'(A)') "========================================="
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
        WRITE(msgBuf,'(A)') "  END SLR_CORR_INIT_FIXED DEBUG "
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
        WRITE(msgBuf,'(A)') "========================================="
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
        endif

      RETURN

      END
CEOC
C----&------------------------------------------------------------------xxxxxxx|
C----&------------------------------------------------------------------xxxxxxx|
#include "SLR_CORR_OPTIONS.h"

CBOP
C     !ROUTINE: SLR_CORR_READPARMS
C     !INTERFACE:
      SUBROUTINE SLR_CORR_READPARMS( myThid )

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE SLR_CORR_READPARMS
C     | o Routine to initialize SLR_CORR variables and constants.
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE
C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "SLR_CORR_PARAMS.h"
#include "SLR_CORR_GRID.h"

C     !INPUT/OUTPUT PARAMETERS:
C     === Routine arguments ===
      INTEGER myThid

#ifdef ALLOW_SLR_CORR

C     !LOCAL VARIABLES:
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      INTEGER iUnit
      INTEGER errCount
      INTEGER i, j
      INTEGER bi, bj, iG, jG, iGm, jGm
      INTEGER OB_indexUnset
CEOP

      NAMELIST /SLR_CORR_PARM01/
     &          useSLR_CORRbalance,
     &          SLR_CORR_balancePeriod,
     &          SLR_target_File

      IF ( .NOT.useSLR_CORR ) THEN
C-    pkg SLR_CORR is not used
        _BEGIN_MASTER(myThid)
C-    Track pkg activation status:
C     print a (weak) warning if data.slr_corr is found
         CALL PACKAGES_UNUSED_MSG( 'useSLR_CORR', ' ', ' ' )
        _END_MASTER(myThid)
        RETURN
      ENDIF

      _BEGIN_MASTER(myThid)

C--   Default flags and values for SLR_CORR
      useSLR_CORRbalance = .TRUE.
      SLR_CORR_balancePeriod = 0. _d 0
      SLR_target_File = ' '

C     Open and read the data.slr_corr file
      WRITE(msgBuf,'(A)') ' SLR_CORR_READPARMS: opening data.slr_corr'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT , myThid )
      CALL OPEN_COPY_DATA_FILE(
     I                          'data.slr_corr', 'SLR_CORR_READPARMS',
     O                          iUnit,
     I                          myThid )

C--   Read parameters from open data file
      READ(UNIT=iUnit,NML=SLR_CORR_PARM01)

      WRITE(msgBuf,'(A)') ' SLR_CORR_READPARMS: finished reading data.slr_corr'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT , myThid )

C--   Close the open data file
#ifdef SINGLE_DISK_IO
      CLOSE(iUnit)
#else
      CLOSE(iUnit,STATUS='DELETE')
#endif /* SINGLE_DISK_IO */

      _END_MASTER(myThid)

C--   Everyone else must wait for the parameters to be loaded
      _BARRIER

#endif /* ALLOW_SLR_CORR */
      RETURN
      END

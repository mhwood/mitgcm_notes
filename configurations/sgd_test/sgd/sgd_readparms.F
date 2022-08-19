#include "SGD_OPTIONS.h"

CBOP
C !ROUTINE: SGD_READPARMS

C !INTERFACE: ==========================================================
      SUBROUTINE SGD_READPARMS( myThid )

C !DESCRIPTION:
C     Initialize SGD parameters, read in data.sgd

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "SGD.h"
#include "PARAMS.h"

C !INPUT PARAMETERS: ===================================================
C  myThid               :: my thread Id number
      INTEGER myThid

C !OUTPUT PARAMETERS: ==================================================
C  none

#ifdef ALLOW_SGD

C !LOCAL VARIABLES: ====================================================
C  iUnit                :: unit number for I/O
C  msgBuf               :: message buffer
      INTEGER iUnit
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      INTEGER, PARAMETER :: debug = 1
      INTEGER m
CEOP

C----&------------------------------------------------------------------xxxxxxx|
C -- Define the variable types in data.sgd
C----&------------------------------------------------------------------xxxxxxx|

      INTEGER sgd_nml_filePrec
      INTEGER sgd_nml_startdates_1(nSGD)
      INTEGER sgd_nml_startdates_2(nSGD)
      _RL sgd_nml_periods(nSGD)
      CHARACTER*30 sgd_nml_mask_fnames(nSGD)
      CHARACTER*30 sgd_nml_volume_flux_fnames(nSGD)


C----&------------------------------------------------------------------xxxxxxx|
C -- Define the namelists in data.sgd
C----&------------------------------------------------------------------xxxxxxx|

      NAMELIST / SGD_PARM01 /
     &     sgd_nml_filePrec

      NAMELIST / SGD_PARM02 /
     &     sgd_nml_mask_fnames,
     &     sgd_nml_volume_flux_fnames,
     &     sgd_nml_startdates_1, sgd_nml_startdates_2,
     &     sgd_nml_periods


C----&------------------------------------------------------------------xxxxxxx|
      if (debug .eq. 1) then
      WRITE(msgBuf,'(A)') "==========================================="
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)') "  BEGIN SGD_READPARMS DEBUG "
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)') "==========================================="
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      endif
C----&------------------------------------------------------------------xxxxxxx|

      IF ( .NOT.useSGD ) THEN
C-    pkg SGD is not used
        _BEGIN_MASTER(myThid)
C     print a (weak) warning if data.sgd is found
         CALL PACKAGES_UNUSED_MSG( 'useSGD', ' ', ' ' )
        _END_MASTER(myThid)
        RETURN
      ENDIF

      _BEGIN_MASTER(myThid)

C----&------------------------------------------------------------------xxxxxxx|
C -- Define the default values for variables in data.sgd
C----&------------------------------------------------------------------xxxxxxx|

C --   Local default flags and values for PRESCRIBE_VEC

      DO m=1,nSGD
            sgd_nml_mask_fnames(m) = ' '
            sgd_nml_volume_flux_fnames(m) = ' '
      ENDDO

      sgd_nml_filePrec  = 32

C----&------------------------------------------------------------------xxxxxxx|

C Open and read the data.sgd file
      WRITE(msgBuf,'(A)') ' SGD_READPARMS: opening data.sgd'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )
      CALL OPEN_COPY_DATA_FILE(
     I                   'data.sgd', 'SGD_READPARMS',
     O                   iUnit,
     I                   myThid )

      WRITE(msgBuf,'(A)') ' SGD_READPARMS: reading SGD_PARM01'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )
      READ(UNIT=iUnit,NML=SGD_PARM01)

      WRITE(msgBuf,'(A)') ' SGD_READPARMS: reading SGD_PARM02'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )
      READ(UNIT=iUnit,NML=SGD_PARM02)

      WRITE(msgBuf,'(A)')
     &  ' SGD_READPARMS: finished reading data.sgd'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT, myThid )

C Close the open data file
#ifdef SINGLE_DISK_IO
      CLOSE(iUnit)
#else
      CLOSE(iUnit,STATUS='DELETE')
#endif /* SINGLE_DISK_IO */

C----&------------------------------------------------------------------xxxxxxx|
C     Fill SGD.h COMMON Block with information provided by data.sgd
C----&------------------------------------------------------------------xxxxxxx|

      DO m=1,nSGD
            sgd_mask_fnames(m) = sgd_nml_mask_fnames(m)
            sgd_volume_flux_fnames(m) = sgd_nml_volume_flux_fnames(m)
            sgd_startdates_1(m) = sgd_nml_startdates_1(m)
            sgd_startdates_2(m) = sgd_nml_startdates_2(m)
            sgd_periods(m) = sgd_nml_periods(m)
      ENDDO

      sgd_filePrec  = sgd_nml_filePrec

C----&------------------------------------------------------------------xxxxxxx|


      _END_MASTER(myThid)
C Everyone else must wait for the parameters to be loaded
      _BARRIER

      if (debug .eq. 1) then
      WRITE(msgBuf,'(A)') "==========================================="
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)') "  END SGD_READPARMS DEBUG "
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)') "==========================================="
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      endif

#endif /* ALLOW_SGD */

      RETURN
      END

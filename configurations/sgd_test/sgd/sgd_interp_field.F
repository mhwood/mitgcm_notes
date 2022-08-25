C------------------------------------------------------------------------------|
C                           sgd_read_field
C------------------------------------------------------------------------------|

#include "SGD_OPTIONS.h"

C----&------------------------------------------------------------------xxxxxxx|
C----&------------------------------------------------------------------xxxxxxx|
CBOP
C !ROUTINE: SGD_READ_FIELD

C !INTERFACE:
      SUBROUTINE SGD_INTERP_FIELD(sgd_id, fac, myThid )

C     !DESCRIPTION:
C     Initialize SGD variables

C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "EESUPPORT.h"
#include "PARAMS.h"
#include "GRID.h"
#include "SGD.h"

C     !INPUT/OUTPUT PARAMETERS:
C     myThid ::  my Thread Id number
      INTEGER sgd_id, i
      INTEGER myThid
      INTEGER mask_len
      _RL fac
      INTEGER, PARAMETER :: debug = 1
      CHARACTER*(MAX_LEN_MBUF) msgBuf
#ifndef ALLOW_USE_MPI
      INTEGER, PARAMETER :: mpiMyId = 0
#endif

CEOP
C------------------------------------------------------------------------------|
CBOC

#ifdef ALLOW_SGD

      if (debug .eq. 1) then
      WRITE(msgBuf,'(A)') "==========================================="
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)') "  BEGIN SGD_INTERP_FIELD DEBUG "
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)') "==========================================="
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      endif

C------------------------------------------------------------------------------|

      mask_len = sgd_numPnts_allproc(sgd_id, mpiMyId+1)
      DO i=1, mask_len
            sgd_fluxes(sgd_id,i) = fac * sgd_fluxes0(sgd_id,i) + 
     &                             (1-fac) * sgd_fluxes1(sgd_id,i)
      ENDDO

C------------------------------------------------------------------------------|

      if (debug .eq. 1) then
      WRITE(msgBuf,'(A)') "==========================================="
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)') "  END SGD_INTERP_FIELD DEBUG "
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)') "==========================================="
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      endif

      RETURN

#endif /* ALLOW_SGD */
      
      END
CEOC
C----&------------------------------------------------------------------xxxxxxx|
C----&------------------------------------------------------------------xxxxxxx|




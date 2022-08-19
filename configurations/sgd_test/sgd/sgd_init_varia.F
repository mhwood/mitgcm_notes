#include "SGD_OPTIONS.h"

CBOP
      SUBROUTINE SGD_INIT_VARIA( myThid )
C     *============================================================*
C     | SUBROUTINE SGD_INIT_VARIA                             |
C     | o Routine to initialize SGD variables.                |
C     *============================================================*
C     | Initialize SGD parameters and variables.              |
C     *============================================================*
      IMPLICIT NONE

C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "SGD.h"

C     === Routine arguments ===
C     myThid -  Number of this instance of SGD_INIT_VARIA
      INTEGER myThid

#ifdef ALLOW_SGD
C     === Local variables ===
C     i, j, k, bi, bj   :: Loop counters
C     errCount          :: error counter
      INTEGER i, j, k, bi, bj, sgd_id
      INTEGER, PARAMETER :: debug = 0
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      LOGICAL read_new_field
      INTEGER iRecord
CEOP
C------------------------------------------------------------------------------|
      if (debug .eq. 1) then
      WRITE(msgBuf,'(A)') "==========================================="
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)') "  BEGIN PRESCRIBE_VEC_READ_FIELDS DEBUG "
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)') "==========================================="
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      endif
C------------------------------------------------------------------------------|

C     Nothing to do here yet

C------------------------------------------------------------------------------|
      if (debug .eq. 1) then
      WRITE(msgBuf,'(A)') "==========================================="
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)') "  END SGD_INIT_VARIA DEBUG "
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)') "==========================================="
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      endif

      RETURN
C------------------------------------------------------------------------------|

#endif /* ALLOW_SGD */

      RETURN
      END

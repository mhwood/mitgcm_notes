#include "SGD_OPTIONS.h"

      SUBROUTINE SGD_CHECK( myThid )
C     /==========================================================\
C     | SUBROUTINE SGD_CHECK                                |
C     | o Validate basic package setup and inter-package         |
C     | dependencies.                                            |
C     \==========================================================/
      IMPLICIT NONE

C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "SGD.h"

C     === Routine arguments ===
C     myThid -  Number of this instance of SGD_CHECK
      INTEGER myThid

#ifdef ALLOW_SGD

C     === Local variables ===
C     msgBuf      - Informational/error meesage buffer
      CHARACTER*(MAX_LEN_MBUF) msgBuf

      _BEGIN_MASTER(myThid)

      WRITE(msgBuf,'(A)') 'SGD_CHECK: #define ALLOW_SGD'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &     SQUEEZE_RIGHT, myThid )

      _END_MASTER(myThid)

#endif /* ALLOW_SGD */

      RETURN
      END

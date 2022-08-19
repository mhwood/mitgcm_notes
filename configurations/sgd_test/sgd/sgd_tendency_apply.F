#include "SGD_OPTIONS.h"
C--  File sgd_tendency_apply.F: Routines to apply SGD tendencies
C--   Contents
C--   o SGD_TENDENCY_APPLY_T
C--   o SGD_TENDENCY_APPLY_S

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP 0
C !ROUTINE: SGD_TENDENCY_APPLY_T

C !INTERFACE:
      SUBROUTINE SGD_TENDENCY_APPLY_T(
     U                    gT_arr,
     I                    k, bi, bj, myTime, myIter, myThid )

C     !DESCRIPTION:
C     Add sgd tendency terms to T (theta) tendency.
C     Routine works for one level at a time.

C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
c#include "GRID.h"
c#include "DYNVARS.h"
#include "SGD.h"

C     !INPUT/OUTPUT PARAMETERS:
C     gT_arr    :: the tendency array
C     k         :: Current vertical level index
C     bi,bj     :: Current tile indices
C     myTime    :: Current time in simulation
C     myIter    :: Current iteration number
C     myThid    :: my Thread Id number
      _RL     gT_arr(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      INTEGER k, bi, bj
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
CEOP

C     !LOCAL VARIABLES:
      INTEGER i, j

C      IF ( applyIcefrontTendT ) THEN
C       DO j=1,sNy
C        DO i=1,sNx
C         IF ( k .LE. k_sgd(i,j,bi,bj) ) THEN  
C          gT_arr(i,j) = gT_arr(i,j)
C     &         + SGD_TendT(i,j,k,bi,bj)
C         ENDIF
C        ENDDO
C       ENDDO
C      ENDIF

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP 0
C !ROUTINE: SGD_TENDENCY_APPLY_S

C !INTERFACE:
      SUBROUTINE SGD_TENDENCY_APPLY_S(
     U                    gS_arr,
     I                    k, bi, bj, myTime, myIter, myThid )

C     !DESCRIPTION:
C     Add sgd tendency terms to S tendency.
C     Routine works for one level at a time.

C     !INPUT PARAMETERS:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
c#include "GRID.h"
c#include "DYNVARS.h"
#include "SGD.h"

C     !INPUT/OUTPUT PARAMETERS:
C     gS_arr    :: the tendency array
C     k         :: Current vertical level index
C     bi,bj     :: Current tile indices
C     myTime    :: Current time in simulation
C     myIter    :: Current iteration number
C     myThid    :: my Thread Id number
      _RL     gS_arr(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      INTEGER k, bi, bj
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
CEOP

C     !LOCAL VARIABLES:
      INTEGER i, j

C      IF ( applyIcefrontTendS ) THEN
C       DO j=1,sNy
C        DO i=1,sNx
C         IF ( k .LE. k_sgd(i,j,bi,bj) ) THEN  
C          gS_arr(i,j) = gS_arr(i,j)
C     &         + SGD_TendS(i,j,k,bi,bj)
C         ENDIF
C        ENDDO
C       ENDDO
C      ENDIF

      RETURN
      END

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
     U                    gT_arr, iMin,iMax,jMin,jMax,
     I                    k, bi, bj, myTime, myIter, myThid )

C     !DESCRIPTION:
C     Add sgd tendency terms to T (theta) tendency.
C     Routine works for one level at a time.

C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "SGD.h"


C     !INPUT/OUTPUT PARAMETERS:
C     gT_arr    :: the tendency array
C     k         :: Current vertical level index
C     bi,bj     :: Current tile indices
C     myTime    :: Current time in simulation
C     myIter    :: Current iteration number
C     myThid    :: my Thread Id number
      _RL     gT_arr(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      INTEGER iMin, iMax, jMin, jMax
      INTEGER k, bi, bj
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      INTEGER, PARAMETER :: debug = 0
CEOP

C     !LOCAL VARIABLES:
      INTEGER sgd_id, s, i, j
      _RL tend_T, cell_vol, prescribed_flux
      _RL tLoc, sLoc, pLoc
#ifndef ALLOW_USE_MPI
      INTEGER, PARAMETER :: mpiMyId = 0
#endif

      _RL SW_TEMP
      EXTERNAL SW_TEMP

C----&------------------------------------------------------------------xxxxxxx|
      if (debug .eq. 1) then
      WRITE(msgBuf,'(A)') "==========================================="
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)') "  BEGIN SGD_TENDENCY_APPLY_T DEBUG "
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)') "==========================================="
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      endif
C----&------------------------------------------------------------------xxxxxxx|
            
C --  If we're going to apply temp
      IF ( sgd_tend_apply_T ) THEN

C --  Loop through each mask 
      DO sgd_id=1,nSGD

C --        For each point in this mask
            DO s=1,sgd_numPnts_allproc(sgd_id, mpiMyId+1)

C --        If this point is in this tile
            IF (k .eq. sgd_point_ij(sgd_id, 3, s)) THEN
            IF (bi .eq. sgd_point_ij(sgd_id, 4, s)) THEN
            IF (bj .eq. sgd_point_ij(sgd_id, 5, s)) THEN

C --        Calculate the tendency of the point

C --        Get the indices      
            i = sgd_point_ij(sgd_id, 1, s)
            j = sgd_point_ij(sgd_id, 2, s)

C --        Get the flux stored from the input file
            prescribed_flux = sgd_fluxes(sgd_id,s)

            if (debug .eq. 1) then
            WRITE(msgBuf,'(A,I5,A,I5,A,I5,A,I5,A,I5)') 
     &             " Applying T tendency at point i=",
     &             i,' j=',j,' k=',k,' bi=',bi,' bj=',bj
            CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
            WRITE(msgBuf,'(A,F8.3,A,F8.3,A)') 
     &             " Flux = ",prescribed_flux,
     &             " m^3/s at temperature = ", sgd_discharge_T," C"
            CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
            endif

C           Make local copies of temperature, salinity and depth (pressure).
            pLoc = ABS(rC(k))
            tLoc = theta(i,j,k,bi,bj)
            sLoc = MAX(salt(i,j,k,bi,bj), 0. _d 0)

C           Turn potential temperature into in-situ temperature [C]
            tLoc = SW_TEMP(sLoc,tLoc,pLoc,0.D0)

C           Calculate the volume of this cell [m^3]
            cell_vol = rA(i,j,bi,bj) * drF(k) * hFacC(i,j,k,bi,bj)

C           Calculate the tendency
            tend_T = ( sgd_discharge_T - tLoc ) * 
     &               prescribed_flux / cell_vol

            WRITE(msgBuf,'(A,F16.12,A,F16.3,A)') 
     &             " dT/dt = ",tend_T,
     &             " C/dt (cell volume = ",
     &             cell_vol," m^3)"
            CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )

C --        Add apply the tendency
            gT_arr(i,j) = gT_arr(i,j) + tend_T

            ENDIF
            ENDIF
            ENDIF

            ENDDO

      ENDDO

      ENDIF

C----&------------------------------------------------------------------xxxxxxx|
      if (debug .eq. 1) then
      WRITE(msgBuf,'(A)') "==========================================="
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)') "  END SGD_TENDENCY_APPLY_T DEBUG "
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)') "==========================================="
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      endif
C----&------------------------------------------------------------------xxxxxxx|

      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP 0
C !ROUTINE: SGD_TENDENCY_APPLY_S

C !INTERFACE:
      SUBROUTINE SGD_TENDENCY_APPLY_S(
     U                    gS_arr, iMin,iMax,jMin,jMax,
     I                    k, bi, bj, myTime, myIter, myThid )

C     !DESCRIPTION:
C     Add sgd tendency terms to S tendency.
C     Routine works for one level at a time.

C     !INPUT PARAMETERS:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "SGD.h"

C     !INPUT/OUTPUT PARAMETERS:
C     gS_arr    :: the tendency array
C     k         :: Current vertical level index
C     bi,bj     :: Current tile indices
C     myTime    :: Current time in simulation
C     myIter    :: Current iteration number
C     myThid    :: my Thread Id number
      _RL     gS_arr(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      INTEGER iMin, iMax, jMin, jMax
      INTEGER k, bi, bj
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      INTEGER, PARAMETER :: debug = 0
CEOP

C     !LOCAL VARIABLES:
      INTEGER sgd_id, s, i, j
      _RL tend_S, cell_vol, prescribed_flux
      _RL sLoc
#ifndef ALLOW_USE_MPI
      INTEGER, PARAMETER :: mpiMyId = 0
#endif

C----&------------------------------------------------------------------xxxxxxx|
      if (debug .eq. 1) then
      WRITE(msgBuf,'(A)') "==========================================="
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)') "  BEGIN SGD_TENDENCY_APPLY_S DEBUG "
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)') "==========================================="
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      endif
C----&------------------------------------------------------------------xxxxxxx|

C --  If we're going to apply temp
      IF ( sgd_tend_apply_T ) THEN

C --  Loop through each mask 
      DO sgd_id=1,nSGD

C --        For each point in this mask
            DO s=1,sgd_numPnts_allproc(sgd_id, mpiMyId+1)

C --        If this point is in this tile
            IF (k .eq. sgd_point_ij(sgd_id, 3, s)) THEN
            IF (bi .eq. sgd_point_ij(sgd_id, 4, s)) THEN
            IF (bj .eq. sgd_point_ij(sgd_id, 5, s)) THEN

C --        Calculate the tendency of the point

C --        Get the indices      
            i = sgd_point_ij(sgd_id, 1, s)
            j = sgd_point_ij(sgd_id, 2, s)

C --        Get the flux stored from the input file
            prescribed_flux = sgd_fluxes(sgd_id,s)

            if (debug .eq. 1) then
            WRITE(msgBuf,'(A,I5,A,I5,A,I5,A,I5,A,I5)') 
     &             " Applying S tendency at point i=",
     &             i,' j=',j,' k=',k,' bi=',bi,' bj=',bj
            CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
            WRITE(msgBuf,'(A,F8.3,A,F8.3,A)') 
     &             " Flux = ",prescribed_flux,
     &             " m^3/s at salinity = ", sgd_discharge_S," kg/m^3"
            CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
            endif

C           Make local copy of salinity 
            sLoc = MAX(salt(i,j,k,bi,bj), 0. _d 0)

C           Calculate the volume of this cell [m^3]
            cell_vol = rA(i,j,bi,bj) * drF(k) * hFacC(i,j,k,bi,bj)
            
C           Calculate the tendency
            tend_S = ( sgd_discharge_S - sLoc ) * 
     &               prescribed_flux / cell_vol

            WRITE(msgBuf,'(A,F16.12,A,F16.3,A)') 
     &             " dS/dt = ",tend_S,
     &             " (kg/m^3)/dt (cell volume = ",
     &             cell_vol," m^3)"
            CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )

C --        Add apply the tendency
            gS_arr(i,j) = gS_arr(i,j) + tend_S

            ENDIF
            ENDIF
            ENDIF

            ENDDO

      ENDDO

      ENDIF

C----&------------------------------------------------------------------xxxxxxx|
      if (debug .eq. 1) then
      WRITE(msgBuf,'(A)') "==========================================="
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)') "  END SGD_TENDENCY_APPLY_S DEBUG "
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)') "==========================================="
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      endif
C----&------------------------------------------------------------------xxxxxxx|

      RETURN
      END

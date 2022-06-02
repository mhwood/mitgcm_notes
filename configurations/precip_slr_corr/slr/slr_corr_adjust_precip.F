#include "SLR_CORR_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C     !ROUTINE: SLR_CORR_ADJUST_PRECIP

C     !INTERFACE:
      SUBROUTINE SLR_CORR_ADJUST_PRECIP( myTime, myIter, myThid )

C     !DESCRIPTION:
C     *==========================================================*
C     | SUBROUTINE SLR_CORR_ADJUST_PRECIP
C     | o Modify precip to ensure match to input slr timeseries
C     *==========================================================*

C     !USES:
      IMPLICIT NONE

C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "SLR_CORR_PARAMS.h"
#include "SLR_CORR_GRID.h"
#ifdef ALLOW_SLR_CORR_SMOOTH_BALANCE
#include "SLR_CORR_BALANCE_SIZE.h"
#include "SLR_CORR_BALANCE_FIELDS.h"
#endif /* ALLOW_SLR_CORR_SMOOTH_BALANCE */
#include "SLR_CORR_FIELDS.h"
#include "FFIELDS.h"


C     !INPUT/OUTPUT PARAMETERS:
      _RL myTime
      INTEGER myIter
      INTEGER myThid
CEOP

#ifdef ALLOW_EXF
#ifdef ALLOW_SLR_CORR_BALANCE

C     !FUNCTIONS:

C     !LOCAL VARIABLES:
C     bi, bj       :: tile indices
C     i,j,k        :: loop indices
C     iB, jB       :: local index of open boundary
C     msgBuf       :: Informational/error message buffer
      INTEGER bi, bj
      INTEGER i, j, k, iB, jB
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      _RL areaOB, areaE, areaW, areaN, areaS, tmpA
      _RL inFlow, flowE, flowW, flowN, flowS
      _RL tileArea(nSx,nSy)
      _RL tileFlow(nSx,nSy)
      _RL tileAreaOB(nSx,nSy)
      _RL tileInFlow(nSx,nSy)
      LOGICAL flag
      _RL netFreshWaterFlux
      _RL shelfIceNetMassFlux
#ifdef ALLOW_SLR_CORR_SMOOTH_BALANCE
      INTEGER n_average_rec
      _RL mean_inFlow
#endif /* ALLOW_SLR_CORR_SMOOTH_BALANCE */      

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

#ifdef ALLOW_DEBUG
      IF (debugMode) CALL DEBUG_ENTER('SLR_CORR_ADJUST_PRECIP',myThid)
#endif

C--   Integrate the transport through each OB
      flag = .FALSE.
      areaOB = 0. _d 0
      inFlow = 0. _d 0
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        tileAreaOB(bi,bj) = 0.
        tileInFlow(bi,bj) = 0.
       ENDDO
      ENDDO

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C--   Calculate a flux throught the surface from the global change of etaH
C     to take into account both freshwater fluxes like precipitation and
C     displacement fluxes from advancing ice sheets or imported sea ice

      netFreshWaterFlux = 0. _d 0
      shelfIceNetMassFlux = 0. _d 0
      IF ( EXFbalanceSurf ) THEN
       DO bj=myByLo(myThid),myByHi(myThid)
        DO bi=myBxLo(myThid),myBxHi(myThid)
         tileFlow(bi,bj) = 0.
         DO j=1,sNy
          DO i=1,sNx
           tileFlow(bi,bj) = tileFlow(bi,bj)
     &          - EmPmR(i,j,bi,bj)
     &          * _rA(i,j,bi,bj) * maskInC(i,j,bi,bj)
          ENDDO
         ENDDO
        ENDDO
       ENDDO
       CALL GLOBAL_SUM_TILE_RL( tileFlow, netFreshWaterFlux, myThid )
       IF ( debugLevel.GE.debLevC ) THEN
        WRITE(msgBuf,'(A,I9,A,1P1E16.8)') 'SLR_CORR_balance (it=',
     &       myIter, ' ) correct for netFreshWaterFlux:',
     &       netFreshWaterFlux
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &       SQUEEZE_RIGHT, myThid )
       ENDIF

      ENDIF

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C--   Calculate a unique velocity correction for all (SLR_CORR_balanceFac>0) OB
C     and correct each OB using corresponding weight factor SLR_CORR_balanceFac

      IF ( flag ) CALL GLOBAL_SUM_TILE_RL( tileAreaOB, areaOB, myThid )
      IF ( areaOB.GT.0. ) THEN
        CALL GLOBAL_SUM_TILE_RL( tileInFlow, inFlow, myThid )
        IF ( debugLevel.GE.debLevC ) THEN
          WRITE(msgBuf,'(A,I9,A,1P2E16.8)') 'SLR_CORR_balance (it=',
     &       myIter, ' ) correct for inFlow:', inFlow, inFlow/areaOB
          CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                        SQUEEZE_RIGHT, myThid )
        ENDIF

        IF ( EXFbalanceSurf ) THEN
         inFlow = inFlow + netFreshWaterFlux*mass2rUnit
#ifdef ALLOW_SHELFICE
         IF ( useShelfIce )
     &        inFlow = inFlow + shelfIceNetMassFlux*mass2rUnit
#endif
        ENDIF
        inFlow = inFlow / areaOB
      ENDIF

C      PRINT *, 'MW: --- CORRECTIONS ----'
C      PRINT *, 'MW: inFlow',inFlow
C      PRINT *, 'MW: areaOB',areaOB
C      PRINT *, 'MW: inFlow/areaOB',inFlow/areaOB

#ifndef ALLOW_SLR_CORR_SMOOTH_BALANCE
      IF ( SLR_CORR_balanceFacE.GE.0. ) flowE =  inFlow*SLR_CORR_balanceFacE
      IF ( SLR_CORR_balanceFacW.GE.0. ) flowW = -inFlow*SLR_CORR_balanceFacW
      IF ( SLR_CORR_balanceFacN.GE.0. ) flowN =  inFlow*SLR_CORR_balanceFacN
      IF ( SLR_CORR_balanceFacS.GE.0. ) flowS = -inFlow*SLR_CORR_balanceFacS
#else
      if (OB_balancePeriod .gt. 0) then
C      PRINT *, 'MW: --- SMOOTH_CORRECTIONS ----'
C        First step is to fill in the flux imbalances vector
C        Before adding a new value, everyone gets shuffled to the right
C        PRINT *, 'SLR_CORR_fluxImbalances',SLR_CORR_fluxImbalances
        DO i=nSLR_CORR_fluxImbalance_rec,2,-1
            SLR_CORR_fluxImbalances(i) = SLR_CORR_fluxImbalances(i-1)
        ENDDO
C        Then, the new flux imbalance is added to the front of the list
        SLR_CORR_fluxImbalances(1) = inFlow
C        And the mean flow is calculated
         n_average_rec = ceiling(OB_balancePeriod/deltaT)
C         PRINT *, 'n_average_rec',n_average_rec
         IF ((myIter - nIter0) .lt. n_average_rec) then
C            PRINT *, 'only using ',myIter - nIter0+1,' values'
C            mean_inFlow = SUM(SLR_CORR_fluxImbalances(1:myIter - nIter0+1))
C            mean_inFlow = mean_inFlow/(myIter - nIter0+1)
             mean_inFlow = 0.
         ELSE
            mean_inFlow = SUM(SLR_CORR_fluxImbalances(1:n_average_rec))
            mean_inFlow = mean_inFlow/n_average_rec
         ENDIF
C         PRINT *, 'SLR_CORR_fluxImbalances',SLR_CORR_fluxImbalances
C        Finally, the mean inFlow is divided among the boundaries

         IF ( SLR_CORR_balanceFacE.GE.0. )
     &       flowE =  mean_inFlow*SLR_CORR_balanceFacE
         IF ( SLR_CORR_balanceFacW.GE.0. )
     &       flowW = -mean_inFlow*SLR_CORR_balanceFacW
         IF ( SLR_CORR_balanceFacN.GE.0. )
     &       flowN =  mean_inFlow*SLR_CORR_balanceFacN
         IF ( SLR_CORR_balanceFacS.GE.0. )
     &       flowS = -mean_inFlow*SLR_CORR_balanceFacS 
       else
         IF ( SLR_CORR_balanceFacE.GE.0. ) flowE =  inFlow*SLR_CORR_balanceFacE
         IF ( SLR_CORR_balanceFacW.GE.0. ) flowW = -inFlow*SLR_CORR_balanceFacW
         IF ( SLR_CORR_balanceFacN.GE.0. ) flowN =  inFlow*SLR_CORR_balanceFacN
         IF ( SLR_CORR_balanceFacS.GE.0. ) flowS = -inFlow*SLR_CORR_balanceFacS 
       endif
#endif /* ALLOW_SLR_CORR_SMOOTH_BALANCE */

C      PRINT *, 'MW: flowE',flowE
C      PRINT *, 'MW: flowW',flowW
C      PRINT *, 'MW: flowS',flowS
C      PRINT *, 'MW: flowN',flowN

      IF ( debugLevel.GE.debLevC  .AND. areaOB.GT.0. ) THEN
          WRITE(msgBuf,'(A,1P2E16.8)')
     &     'SLR_CORR_balance correction to OBEu,OBWu:', flowE, flowW
          CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                        SQUEEZE_RIGHT, myThid )
          WRITE(msgBuf,'(A,1P2E16.8)')
     &     'SLR_CORR_balance correction to OBNv,OBSv:', flowN, flowS
          CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                        SQUEEZE_RIGHT, myThid )
      ENDIF

c     IF ( .NOT.useEXFbalance ) RETURN

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C--   Add correction:

#ifdef ALLOW_SLR_CORR_EAST
      IF ( SLR_CORR_balanceFacE.NE.0. ) THEN
       DO bj=myByLo(myThid),myByHi(myThid)
        DO bi=myBxLo(myThid),myBxHi(myThid)
         IF ( tileHasOBE(bi,bj) ) THEN
          DO k=1,Nr
           DO j=1-OLy,sNy+OLy
            iB = OB_Ie(j,bi,bj)
            IF ( iB.NE.OB_indexNone ) THEN
             OBEu(j,k,bi,bj) = OBEu(j,k,bi,bj)
     &                       + flowE*maskW(iB,j,k,bi,bj)
            ENDIF
           ENDDO
          ENDDO
         ENDIF
        ENDDO
       ENDDO
      ENDIF
#endif /* ALLOW_SLR_CORR_EAST */

#ifdef ALLOW_SLR_CORR_WEST
      IF ( SLR_CORR_balanceFacW.NE.0. ) THEN
       DO bj=myByLo(myThid),myByHi(myThid)
        DO bi=myBxLo(myThid),myBxHi(myThid)
         IF ( tileHasOBW(bi,bj) ) THEN
          DO k=1,Nr
           DO j=1-OLy,sNy+OLy
            iB = OB_Iw(j,bi,bj)
            IF ( iB.NE.OB_indexNone ) THEN
             OBWu(j,k,bi,bj) = OBWu(j,k,bi,bj)
     &                       + flowW*maskW(1+iB,j,k,bi,bj)
            ENDIF
           ENDDO
          ENDDO
         ENDIF
        ENDDO
       ENDDO
      ENDIF
#endif /* ALLOW_SLR_CORR_WEST */

#ifdef ALLOW_SLR_CORR_NORTH
      IF ( SLR_CORR_balanceFacN.NE.0. ) THEN
       DO bj=myByLo(myThid),myByHi(myThid)
        DO bi=myBxLo(myThid),myBxHi(myThid)
         IF ( tileHasOBN(bi,bj) ) THEN
          DO k=1,Nr
           DO i=1-OLx,sNx+OLx
            jB = OB_Jn(i,bi,bj)
            IF ( jB.NE.OB_indexNone ) THEN
             OBNv(i,k,bi,bj) = OBNv(i,k,bi,bj)
     &                       + flowN*maskS(i,jB,k,bi,bj)
            ENDIF
           ENDDO
          ENDDO
         ENDIF
        ENDDO
       ENDDO
      ENDIF
#endif /* ALLOW_SLR_CORR_NORTH */

#ifdef ALLOW_SLR_CORR_SOUTH
      IF ( SLR_CORR_balanceFacS.NE.0. ) THEN
       DO bj=myByLo(myThid),myByHi(myThid)
        DO bi=myBxLo(myThid),myBxHi(myThid)
         IF ( tileHasOBS(bi,bj) ) THEN
          DO k=1,Nr
           DO i=1-OLx,sNx+OLx
            jB = OB_Js(i,bi,bj)
            IF ( jB.NE.OB_indexNone ) THEN
             OBSv(i,k,bi,bj) = OBSv(i,k,bi,bj)
     &                       + flowS*maskS(i,1+jB,k,bi,bj)
            ENDIF
           ENDDO
          ENDDO
         ENDIF
        ENDDO
       ENDDO
      ENDIF
#endif /* ALLOW_SLR_CORR_SOUTH */

#ifdef ALLOW_DEBUG
      IF (debugMode) CALL DEBUG_LEAVE('SLR_CORR_ADJUST_PRECIP',myThid)
#endif

#endif /* ALLOW_SLR_CORR_BALANCE */
#endif /* ALLOW_EXF */

      RETURN
      END

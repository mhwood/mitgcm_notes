#include "SGD_OPTIONS.h"

CBOP
      SUBROUTINE SGD_INIT_FIXED( myThid )
C     *============================================================*
C     | SUBROUTINE SGD_INIT_FIXED
C     | o Routine to initialize SGD parameters and variables.
C     *============================================================*
C     | Initialize SGD parameters and variables.
C     *============================================================*
      IMPLICIT NONE

C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "SGD.h"

C     === Routine arguments ===
C     myThid        :: my Thread Id Number
      INTEGER myThid

#ifdef ALLOW_SGD
C     === Local variables ===
C     i, j, k, bi, bj   :: Loop counters
C     errCount          :: error counter
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      INTEGER sgd_id, i, j, k, bi, bj, p
      INTEGER, PARAMETER :: debug = 1
CEOP

      INTEGER mask_name_len
      CHARACTER*50 tmpstr
      CHARACTER*50 mask_filename
      _RL mask_buffer(nSGD,1-Olx:sNx+Olx,1-Oly:sNy+Oly,1:Nr,nSx,nSy)
      _RL mask_start_time


#ifdef ALLOW_CAL
      INTEGER date_array(4), difftime(4), yearStartDate(4)
      INTEGER gcm_startdate(4)
#endif /* ALLOW_CAL */

      if (debug .eq. 1) then
      WRITE(msgBuf,'(A)') "==========================================="
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)') "  BEGIN SGD_INIT_FIXED DEBUG "
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)') "==========================================="
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      endif

C----&------------------------------------------------------------------xxxxxxx|

C----&------------------------------------------------------------------xxxxxxx|
C     Fill SGD.h COMMON Block with default 
C          values for the vector variables
C----&------------------------------------------------------------------xxxxxxx|


      DO sgd_id=1, nSGD

            sgd_nPoints(sgd_id) = 0

            DO p=1,nPx*nPy
                  sgd_numPnts_allproc(sgd_id, p) = 0
            ENDDO

            DO bj = myByLo(myThid), myByHi(myThid)
            DO bi = myBxLo(myThid), myBxHi(myThid)
                  DO j=1-OLy,sNy+OLy
                  DO i=1-OLx,sNx+OLx
                        DO k=1,Nr
                        sgd_subMask(sgd_id,i,j,k,bi,bj) = 0.0
                        ENDDO
                  ENDDO
                  ENDDO
            ENDDO
            ENDDO

            DO i=1,sNx*sNy
            DO p=1,nPx*nPy
                  sgd_mask_index_list(sgd_id, p, i) = 0
            ENDDO
            ENDDO

            DO i=1,sNx*sNy
                  sgd_point_ij(sgd_id, 1, i) = 0
                  sgd_point_ij(sgd_id, 2, i) = 0
                  sgd_point_ij(sgd_id, 3, i) = 0
                  sgd_point_ij(sgd_id, 4, i) = 0
                  sgd_point_ij(sgd_id, 5, i) = 0
            ENDDO

            DO i=1,max_SGD_points
                  sgd_fluxes(sgd_id,i) = 0.0
                  sgd_fluxes0(sgd_id,i) = 0.0
                  sgd_fluxes1(sgd_id,i) = 0.0
            ENDDO

      ENDDO

#ifdef ALLOW_CAL
      DO i=1,nSGD
      CALL CAL_FULLDATE( sgd_startdates_1(i), sgd_startdates_2(i),
     &                   date_array, myThid )
          CALL cal_getdate( nIter0, startTime, gcm_startdate, myThid )
          CALL cal_TimePassed( gcm_startdate, date_array, difftime,
     I                         myThid )
          CALL cal_ToSeconds ( difftime, mask_start_time, myThid )
          mask_start_time = startTime  + mask_start_time
          sgd_starttimes(i) = mask_start_time
      ENDDO
#endif /* ALLOW_CAL */

      PRINT *,sgd_starttimes

C----&------------------------------------------------------------------xxxxxxx|
C     Read in the masks
C----&------------------------------------------------------------------xxxxxxx|

C     Looping through each vector mask
      DO sgd_id=1, nSGD

      mask_filename = sgd_mask_fnames(sgd_id)

      if (debug .eq. 1) then
        
        WRITE(msgBuf,'(A,A,A)')  'Identifying sgd mask points for ',
     &  ' mask ', mask_filename
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      endif

C     Counting number of valid characters to trim trailing whitespaces
      mask_name_len = 0
      tmpstr = sgd_mask_fnames(sgd_id)
      DO i=1, len(tmpstr)
        IF (tmpstr(i:i) .ne. ' ') then
            mask_name_len = mask_name_len + 1
        ENDIF
      ENDDO
      WRITE(mask_filename,'(A)') 
     &     tmpstr(1:mask_name_len)

C       Read mask files and assign portion of mask to each process
        CALL READ_REC_XYZ_RL( mask_filename,
     &               mask_buffer(sgd_id,:,:,:,:,:), 1, 0, myThid )

C       Filling in halo regions of sub masks
        _EXCH_XY_RL(mask_buffer(sgd_id,:,:,:,:,:), myThid)

        DO bj = myByLo(myThid), myByHi(myThid)
            DO bi = myBxLo(myThid), myBxHi(myThid)
                  DO j=1-OLy,sNy+OLy
                  DO i=1-OLx,sNx+OLx
                        DO k=1,Nr
                        sgd_subMask(sgd_id,i,j,k,bi,bj)
     &                      = mask_buffer(sgd_id,i,j,k,bi,bj)
                        ENDDO
                  ENDDO
                  ENDDO
            ENDDO
        ENDDO

C       Mask a list of where the mask points are located
        CALL IDENTIFY_SGD_POINTS(sgd_subMask,sgd_mask_index_list,
     &                              sgd_id, myThid)

      ENDDO

C----&------------------------------------------------------------------xxxxxxx|


      if (debug .eq. 1) then
      WRITE(msgBuf,'(A)') "==========================================="
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)') "  END SGD_INIT_FIXED DEBUG "
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)') "==========================================="
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      endif

#endif /* ALLOW_SGD */

      RETURN
      END



#ifdef ALLOW_SGD
C----&------------------------------------------------------------------xxxxxxx|
C----&------------------------------------------------------------------xxxxxxx|
CBOP
C !ROUTINE: IDENTIFY_SGD_POINTS

C !INTERFACE:
      SUBROUTINE IDENTIFY_SGD_POINTS(mask_buffer,index_buffer,
     &       sgd_id, myThid)
 
C     !DESCRIPTION:
C     Locate vector mask points within each subtile
       IMPLICIT NONE
 
#include "SIZE.h"
#include "EEPARAMS.h"
#include "EESUPPORT.h"
#include "PARAMS.h"
#include "SGD.h"
 
C     !INPUT PARAMETERS:
      _RL mask_buffer(nSGD,1-Olx:sNx+Olx,1-Oly:sNy+Oly,1:Nr,nSx,nSy)
      INTEGER index_buffer(nSGD, nPx*nPy, sNy*sNx*Nr)
      INTEGER sgd_id, myThid
 
C     !LOCAL VARIABLES:     
      INTEGER counter, k, i, j, bi, bj
      INTEGER mask_val
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      INTEGER, PARAMETER :: debug = 1
#ifndef ALLOW_USE_MPI
      INTEGER, PARAMETER :: mpiMyId = 0
#endif
 
CEOP
C----&------------------------------------------------------------------xxxxxxx|
CBOC
 
      if (debug .eq. 1) then
      WRITE(msgBuf,'(A)') "==========================================="
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)') "  BEGIN IDENTIFY_SGD_POINTS DEBUG "
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)') "==========================================="
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      endif

      counter = 1
      DO bj = myByLo(myThid), myByHi(myThid)
      DO bi = myBxLo(myThid), myBxHi(myThid)
      DO k=1,Nr
      DO j=1,sNy
      DO i=1,sNx
      
C            if (k.eq.48) then
C            PRINT*,sgd_id,k,i,j,bi,bj,mask_buffer(sgd_id,k,i,j,bi,bj)
C            endif
            
            if (mask_buffer(sgd_id,i,j,k,bi,bj) .ge. 1.0 ) then

            mask_val = INT(mask_buffer(sgd_id,i,j,k,bi,bj))

C           Save the mask index for this processor
            sgd_point_ij(sgd_id,1,counter) = i
            sgd_point_ij(sgd_id,2,counter) = j
            sgd_point_ij(sgd_id,3,counter) = k
            sgd_point_ij(sgd_id,4,counter) = bi
            sgd_point_ij(sgd_id,5,counter) = bj

            index_buffer(sgd_id,mpiMyId+1,counter) = mask_val

            counter = counter + 1
            end if
            ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
 
      sgd_numPnts_allproc(sgd_id,mpiMyId+1) = counter-1

      WRITE(msgBuf,'(A,I5,A,I5,A,I5)') 
     &   "sgd_numPnts_allproc(",sgd_id,",",mpiMyId+1,
     &   ") = ",sgd_numPnts_allproc(sgd_id,mpiMyId+1)
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &        SQUEEZE_RIGHT, myThid )

      if ((debug .eq. 1)) then
      DO i=1,sgd_numPnts_allproc(sgd_id,mpiMyId+1)
            WRITE(msgBuf,'(A,I5,A,I5,A,I5,A,I5,A,I5)') 
     &     "sgd_point_ij(",sgd_id,",:,",i,
     &     ") = ",sgd_point_ij(sgd_id,1,i),',',
     &     sgd_point_ij(sgd_id,2,i), ',',
     &     sgd_point_ij(sgd_id,3,i)
            CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &               SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A,I5,A,I5,A,I5,A,I5)') 
     &     "mask_index(",sgd_id,",",mpiMyId+1,',',i,
     &     ") = ",index_buffer(sgd_id,mpiMyId+1,i)
            CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &               SQUEEZE_RIGHT, myThid )
      ENDDO
      endif


      if (debug .eq. 1) then
      WRITE(msgBuf,'(A)') "==========================================="
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)') "  END IDENTIFY_SGD_POINTS DEBUG "
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)') "==========================================="
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      endif
 
      END

CEOC
C----&------------------------------------------------------------------xxxxxxx|
C----&------------------------------------------------------------------xxxxxxx|
#endif /* ALLOW_SGD */


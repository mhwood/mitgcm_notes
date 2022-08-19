C------------------------------------------------------------------------------|
C                           sgd_read_field
C------------------------------------------------------------------------------|

#include "SGD_OPTIONS.h"

C----&------------------------------------------------------------------xxxxxxx|
C----&------------------------------------------------------------------xxxxxxx|
CBOP
C !ROUTINE: SGD_READ_FIELD

C !INTERFACE:
      SUBROUTINE SGD_READ_FIELD(file_name, sgd_id,
     I                          iRecord, myThid )

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
      CHARACTER(*) file_name
      INTEGER sgd_id
      INTEGER myThid
      INTEGER iRecord


C !FUNCTIONS
      INTEGER  MDS_RECLEN
      EXTERNAL MDS_RECLEN

C     !LOCAL VARIABLES:
C      INTEGER i,j,k,l, bi,bj
      INTEGER, PARAMETER :: debug = 1
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      INTEGER xSize, ySize
      INTEGER length_of_rec,rcl
      INTEGER global_mask_len
      INTEGER irec, dUnit, i, k
      Real*4 global_buff_r4(max_SGD_points)
      Real*8 global_buff_r8(max_SGD_points)
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
      WRITE(msgBuf,'(A)') "  BEGIN SGD_READ_FIELD DEBUG "
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)') "==========================================="
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      endif

C------------------------------------------------------------------------------|

C - fill in the buffers with zeros
      IF (sgd_filePrec .EQ. 32) THEN
        DO i=1,max_SGD_points
            global_buff_r4(i) = 0.0
        ENDDO
      ELSE
        DO i=1,max_SGD_points
            global_buff_r8(i) = 0.0
        ENDDO
      ENDIF

      global_mask_len = 0
      DO i=1,nPx*nPy
            global_mask_len = global_mask_len + 
     &       sgd_numPnts_allproc(nSGD, nPx*nPy)
      ENDDO

      if (debug .eq. 1) then
      WRITE(msgBuf,'(A,A)') "  Reading from ",file_name
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A,I5,A,I5,A,I5,A,I5)') "    mask_len = ",
     &                   sgd_numPnts_allproc(sgd_id, mpiMyId+1),
     &                  ",      global_mask_len =",global_mask_len,
     &                            ",    iRecord =",iRecord  
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      endif

C  - This part reads in the GLOBAL variable for this timestep (iRecord number)

      CALL MDSFINDUNIT( dUnit, myThid )
      length_of_rec = MDS_RECLEN(sgd_filePrec, global_mask_len, myThid)
      
      OPEN( dUnit, file=file_name, status='old',
     &     access='direct', recl=length_of_rec )
      irec = iRecord
      IF (sgd_filePrec .EQ. precFloat32) THEN
      READ(dUnit,rec=irec) (global_buff_r4(i), i=1,global_mask_len)
      ELSE
      READ(dUnit,rec=irec) (global_buff_r8(i), i=1,global_mask_len)
      ENDIF
      CLOSE( dUnit )

C  - This part then subsets the global variable onto the part for this tile
C     
C      IF (sgd_filePrec .EQ. precFloat32) THEN
C            CALL PASS_GLOBAL_R4_to_RL_SUBSET(global_buff_r4,
C     &                                       sgd_id, fld,
C     &                                       mask_len, n_levels)
C      ELSE
C            CALL PASS_GLOBAL_R8_to_RL_SUBSET(global_buff_r8,
C     &                                       sgd_id, fld,
C     &                                       mask_len, n_levels)
C      ENDIF

C------------------------------------------------------------------------------|

      if (debug .eq. 1) then
      WRITE(msgBuf,'(A)') "==========================================="
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)') "  END SGD_READ_FIELD DEBUG "
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



! C----&------------------------------------------------------------------xxxxxxx|
! C----&------------------------------------------------------------------xxxxxxx|
!        SUBROUTINE PASS_GLOBAL_R4_to_RL_SUBSET(global_buff_r4,
!      &                                        sgd_id, fld,
!      &                                        mask_len, n_levels)

!       IMPLICIT NONE
! #include "SIZE.h"
! #include "EEPARAMS.h"
! #include "EESUPPORT.h"
! #include "SGD.h"

!          Real*4 global_buff_r4(max_SGD_points)
!          INTEGER sgd_id, fld, mask_len, n_levels
!          INTEGER global_index
!          INTEGER i, k
!          INTEGER, PARAMETER :: debug = 1
! #ifndef ALLOW_USE_MPI
!          INTEGER, PARAMETER :: mpiMyId = 0
! #endif
! C----&------------------------------------------------------------------xxxxxxx|
! C----&------------------------------------------------------------------xxxxxxx|
           
!            DO i=1, mask_len
!             global_index = sgd_mask_index_list(sgd_id,mpiMyId+1,i)
!             DO k=1, n_levels
!             sgd_subFields(sgd_id,fld,i,k) = 
!      &       global_buff_r4(global_index,k)
!             ENDDO
!            ENDDO

!        END
! C----&------------------------------------------------------------------xxxxxxx|
! C----&------------------------------------------------------------------xxxxxxx|


! C----&------------------------------------------------------------------xxxxxxx|
! C----&------------------------------------------------------------------xxxxxxx|
!        SUBROUTINE PASS_GLOBAL_R8_to_RL_SUBSET(global_buff_r8,
!      &                                        sgd_id, fld,
!      &                                        mask_len, n_levels)

!       IMPLICIT NONE

! #include "SIZE.h"
! #include "EEPARAMS.h"
! #include "EESUPPORT.h"
! #include "SGD.h"
 
!           Real*8 global_buff_r8(sNx*sNy*nPx*nPy,Nr)
!           INTEGER sgd_id, fld, mask_len, n_levels
!           INTEGER global_index
!           INTEGER i, k
!           INTEGER, PARAMETER :: debug = 0
! #ifndef ALLOW_USE_MPI
!           INTEGER, PARAMETER :: mpiMyId = 0
!  #endif
! C----&------------------------------------------------------------------xxxxxxx|
! C----&------------------------------------------------------------------xxxxxxx|
!       DO i=1, mask_len
!       global_index = sgd_mask_index_list(sgd_id,mpiMyId+1,i)
!       DO k=1, n_levels
!             sgd_fluxes(sgd_id,i) = 0.0 = 
!      &       global_buff_r8(global_index,k)
!       ENDDO
!       ENDDO
 
!         END
! C----&------------------------------------------------------------------xxxxxxx|
! C----&------------------------------------------------------------------xxxxxxx|
 





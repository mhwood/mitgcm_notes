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
     I                          iRecord, read_array, myThid )

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
      INTEGER read_array
      INTEGER length_of_rec
      INTEGER global_mask_len, mask_len
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
      mask_len = sgd_numPnts_allproc(sgd_id, mpiMyId+1)
     
      IF (sgd_filePrec .EQ. precFloat32) THEN
            CALL PASS_GLOBAL_R4_to_RL_SUBSET(global_buff_r4,
     &                                       sgd_id, read_array, 
     &                                       mask_len, myThid)
      ELSE
            CALL PASS_GLOBAL_R8_to_RL_SUBSET(global_buff_r8,
     &                                       sgd_id, read_array, 
     &                                       mask_len, myThid)
      ENDIF

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



C----&------------------------------------------------------------------xxxxxxx|
C----&------------------------------------------------------------------xxxxxxx|
        SUBROUTINE PASS_GLOBAL_R4_to_RL_SUBSET(global_buff_r4,
     &                                         sgd_id, read_array, 
     &                                         mask_len, myThid)
       IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "EESUPPORT.h"
#include "SGD.h"

          Real*4 global_buff_r4(max_SGD_points)
          CHARACTER*(MAX_LEN_MBUF) msgBuf
          INTEGER myThid
          INTEGER sgd_id, fld, mask_len, n_levels
          INTEGER global_index, read_array
          INTEGER i, k
          INTEGER, PARAMETER :: debug = 0
#ifndef ALLOW_USE_MPI
            INTEGER, PARAMETER :: mpiMyId = 0
#endif
C----&------------------------------------------------------------------xxxxxxx|
C----&------------------------------------------------------------------xxxxxxx|
           
      if (debug .eq. 1) then
      WRITE(msgBuf,'(A)') "==========================================="
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)') "  BEGIN PASS_GLOBAL_R4_to_RL_SUBSET DEBUG "
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)') "==========================================="
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      endif

            DO i=1, mask_len
             global_index = sgd_mask_index_list(sgd_id,mpiMyId+1,i)
             if (read_array.eq.0) then
             sgd_fluxes0(sgd_id,i) = global_buff_r4(global_index)
             endif
             if (read_array.eq.1) then
             sgd_fluxes1(sgd_id,i) = global_buff_r4(global_index)
             endif
           ENDDO

      if (debug .eq. 1) then
      WRITE(msgBuf,'(A)') "==========================================="
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)') "  END PASS_GLOBAL_R4_to_RL_SUBSET DEBUG "
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)') "==========================================="
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      endif

        END
C----&------------------------------------------------------------------xxxxxxx|
C----&------------------------------------------------------------------xxxxxxx|


C----&------------------------------------------------------------------xxxxxxx|
C----&------------------------------------------------------------------xxxxxxx|
        SUBROUTINE PASS_GLOBAL_R8_to_RL_SUBSET(global_buff_r8,
     &                                        sgd_id, read_array, 
     &                                        mask_len, myThid)
       IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "EESUPPORT.h"
#include "SGD.h"

          Real*8 global_buff_r8(max_SGD_points)
          CHARACTER*(MAX_LEN_MBUF) msgBuf
          INTEGER myThid
          INTEGER sgd_id, fld, mask_len, n_levels
          INTEGER global_index, read_array
          INTEGER i, k
          INTEGER, PARAMETER :: debug = 1
#ifndef ALLOW_USE_MPI
            INTEGER, PARAMETER :: mpiMyId = 0
#endif
C----&------------------------------------------------------------------xxxxxxx|
C----&------------------------------------------------------------------xxxxxxx|
           
            DO i=1, mask_len
             global_index = sgd_mask_index_list(sgd_id,mpiMyId+1,i)
             if (read_array.eq.0) then
             sgd_fluxes0(sgd_id,i) = global_buff_r8(global_index)
             endif
             if (read_array.eq.1) then
             sgd_fluxes1(sgd_id,i) = global_buff_r8(global_index)
             endif
           ENDDO

        END
C----&------------------------------------------------------------------xxxxxxx|
C----&------------------------------------------------------------------xxxxxxx|
 





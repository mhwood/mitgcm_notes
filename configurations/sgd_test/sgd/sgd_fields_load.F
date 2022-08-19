#include "SGD_OPTIONS.h"

CBOP
      SUBROUTINE SGD_FIELDS_LOAD( myTime, myIter, myThid )
C     *============================================================*
C     | SGD_FIELDS_LOAD                                              |
C     | o Routine to load in SGD fluxes.                           |
C     *============================================================*
      IMPLICIT NONE

C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "SGD.h"

C     === Routine arguments ===
C     myThid -  Number of this instance of SGD_FIELDS_LOAD
      _RL     myTime
      INTEGER myIter
      INTEGER myThid

#ifdef ALLOW_SGD
C     === Local variables ===
C     i, j, k, bi, bj   :: Loop counters
C     errCount          :: error counter
      INTEGER i, j, k, bi, bj, sgd_id
      INTEGER, PARAMETER :: debug = 1
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      LOGICAL read_new_field
      INTEGER iRecord
      INTEGER flux_name_len
      CHARACTER*50 flux_filename

      _RL fldStartTime, fldPeriod
      LOGICAL first, changed
      INTEGER count0, count1, year0, year1
      _RL     fac
      _RL, PARAMETER :: fldRepeatCycleDummy = 0.0
      LOGICAL, PARAMETER :: usefldyearlyfieldsDummy = .FALSE.

CEOP
C------------------------------------------------------------------------------|
      if (debug .eq. 1) then
      WRITE(msgBuf,'(A)') "==========================================="
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)') "  BEGIN SGD_FIELDS_LOAD DEBUG "
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)') "==========================================="
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      endif
C------------------------------------------------------------------------------|

C     Loop through each mask to read in the fluxes (if necessary)
      DO sgd_id=1,nSGD

C           Determine whether a new field should be read from the file
C           based on the date
            read_new_field = .TRUE.
            
            IF (read_new_field) THEN

C           Counting number of valid characters to trim trailing whitespaces 
C           around the flux file name
            flux_name_len = 0
            DO i=1, len(sgd_volume_flux_fnames(sgd_id))
                  IF (sgd_volume_flux_fnames(sgd_id)(i:i) .ne. ' ') then
                  flux_name_len = flux_name_len + 1
                  ENDIF
            ENDDO
            flux_filename = sgd_volume_flux_fnames(sgd_id)

            if (debug .eq. 1) then
            WRITE(msgBuf,'(A,A,A)')  'Reading new fluxes from  ',
     &                               ' file ', flux_filename
            CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
            endif

            print *,'sgd_periods(sgd_id)',sgd_periods(sgd_id)

C           Determine the record to read from the file
            CALL EXF_GetFFieldRec(
     I               sgd_starttimes(sgd_id), sgd_periods(sgd_id), 
     I               fldRepeatCycleDummy,
     I               'SGD', usefldyearlyfieldsDummy,
     O               fac, first, changed,
     O               count0, count1, year0, year1,
     I               myTime, myIter, myThid )
            iRecord = myIter+1

            PRINT *, 'fac, first, changed', fac, first, changed
            PRINT *, 'count0, count1, year0, year1',
     &        count0, count1, year0, year1

            CALL SGD_READ_FIELD(flux_filename, sgd_id, iRecord, myThid )

            ENDIF

      ENDDO

C------------------------------------------------------------------------------|
      if (debug .eq. 1) then
      WRITE(msgBuf,'(A)') "==========================================="
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      WRITE(msgBuf,'(A)') "  END SGD_FIELDS_LOAD DEBUG "
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

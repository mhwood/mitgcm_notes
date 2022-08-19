#ifdef ALLOW_SGD

C------------------------------------------------------------------------------|
C                           SGD.h
C------------------------------------------------------------------------------|

#include "SGD_SIZE.h"

CBOP
C !ROUTINE: SGD.h

C !DESCRIPTION: \bv
C     *==========================================================*
C     | SGD.h                                                    |
C     | o Basic header subglacial discharge package.             |
C     |   Contains all SGD field declarations.                   |
C     *==========================================================*
C \ev

C-----------------------------------------------------------------------
C--   Constants that can be set in data.sgd
C-----------------------------------------------------------------------
CEOP


C     These are the file names of the input masks
C     These masks identify where the sgd outlet locatations are
      CHARACTER*30 sgd_mask_fnames(nSGD)

C     These are the file names of the volume flux values
      CHARACTER*30 sgd_volume_flux_fnames(nSGD)

C     This is the number of points in each mask
      INTEGER sgd_nPoints(nSGD)

C     These store a list of mask points that each proc has
      INTEGER sgd_numPnts_allproc(nSGD, nPx*nPy)

C     This is where the input masks are stored after they are read in
      _RL sgd_subMask(nSGD,1-Olx:sNx+Olx,1-Oly:sNy+Oly,1:Nr,nSx,nSy)

C     These are dictionaries linking the mask counter to its proc
C           tile, and the row/col within that tile
      INTEGER sgd_mask_index_list(nSGD, nPx*nPy, sNy*sNx*Nr)

C     This is where the local depth row and col of the location in the tile are stored
      INTEGER sgd_point_ij(nSGD, 5, sNy*sNx)

C     These are some i/o parameters to help read in the fields
      INTEGER sgd_filePrec
      INTEGER sgd_startdates_1(nSGD)
      INTEGER sgd_startdates_2(nSGD)
      _RL     sgd_periods(nSGD)
      _RL     sgd_starttimes(nSGD)

C     This is where the input masks are stored after they are read in
      _RL sgd_fluxes(nSGD,max_SGD_points)
      _RL sgd_fluxes0(nSGD,max_SGD_points)
      _RL sgd_fluxes1(nSGD,max_SGD_points)

C-----------------------------------------------------------------------

      COMMON / DIAG_VEC_VARS_I /
     &     sgd_nPoints, sgd_numPnts_allproc, 
     &     sgd_mask_index_list, sgd_point_ij,
     &     sgd_startdates_1, sgd_startdates_2,
     &     sgd_filePrec

      COMMON / DIAG_VEC_VARS_R /
     &     sgd_subMask, sgd_fluxes,
     &     sgd_fluxes0, sgd_fluxes1,
     &     sgd_periods, sgd_starttimes

      COMMON / DIAG_VEC_VARS_C /
     &     sgd_mask_fnames,
     &     sgd_volume_flux_fnames

#endif /* ALLOW_SGD */

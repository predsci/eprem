/* The Earth-Moon-Mars Radiation Environment Module (EMMREM) software is */
/* free software; you can redistribute and/or modify the EMMREM sotware */
/* or any part of the EMMREM software under the terms of the GNU General */
/* Public License (GPL) as published by the Free Software Foundation; */
/* either version 2 of the License, or (at your option) any later */
/* version. Software that uses any portion of the EMMREM software must */
/* also be released under the GNU GPL license (version 2 of the GNU GPL */
/* license or a later version). A copy of this GNU General Public License */
/* may be obtained by writing to the Free Software Foundation, Inc., 59 */
/* Temple Place, Suite 330, Boston MA 02111-1307 USA or by viewing the */
/* license online at http://www.gnu.org/copyleft/gpl.html. */

#ifndef READMAS_H
#define READMAS_H
#include <hdf5.h>
#include <hdf5_hl.h>
#include "cubeShellStruct.h"

extern Scalar_t *masTime;
extern Scalar_t *masHelTime;

extern Scalar_t phiOffset;
extern Scalar_t phiHelOffset;

extern Index_t masFileIndex0;
extern Index_t masFileIndex1;
extern Index_t masHelFileIndex0;
extern Index_t masHelFileIndex1;

extern Index_t masMallocFlag;
extern Index_t masEqFileFlag;
extern Index_t masResolutionChangeFlag;

extern int32_t masDimMin_0[1];
extern int32_t masDimMin_1[1];

extern int32_t masBppDimMax_0[1];
extern int32_t masBptDimMax_0[1];
extern int32_t masBprDimMax_0[1];
extern int32_t masHelBppDimMax_0[1];
extern int32_t masHelBptDimMax_0[1];
extern int32_t masHelBprDimMax_0[1];

extern int32_t masBppDimMax_1[1];
extern int32_t masBptDimMax_1[1];
extern int32_t masBprDimMax_1[1];
extern int32_t masHelBppDimMax_1[1];
extern int32_t masHelBptDimMax_1[1];
extern int32_t masHelBprDimMax_1[1];

extern int32_t masBtpDimMax_0[1];
extern int32_t masBttDimMax_0[1];
extern int32_t masBtrDimMax_0[1];
extern int32_t masHelBtpDimMax_0[1];
extern int32_t masHelBttDimMax_0[1];
extern int32_t masHelBtrDimMax_0[1];

extern int32_t masBtpDimMax_1[1];
extern int32_t masBttDimMax_1[1];
extern int32_t masBtrDimMax_1[1];
extern int32_t masHelBtpDimMax_1[1];
extern int32_t masHelBttDimMax_1[1];
extern int32_t masHelBtrDimMax_1[1];

extern int32_t masBrpDimMax_0[1];
extern int32_t masBrtDimMax_0[1];
extern int32_t masBrrDimMax_0[1];
extern int32_t masHelBrpDimMax_0[1];
extern int32_t masHelBrtDimMax_0[1];
extern int32_t masHelBrrDimMax_0[1];

extern int32_t masBrpDimMax_1[1];
extern int32_t masBrtDimMax_1[1];
extern int32_t masBrrDimMax_1[1];
extern int32_t masHelBrpDimMax_1[1];
extern int32_t masHelBrtDimMax_1[1];
extern int32_t masHelBrrDimMax_1[1];

extern int32_t masVppDimMax_0[1];
extern int32_t masVptDimMax_0[1];
extern int32_t masVprDimMax_0[1];
extern int32_t masHelVppDimMax_0[1];
extern int32_t masHelVptDimMax_0[1];
extern int32_t masHelVprDimMax_0[1];

extern int32_t masVppDimMax_1[1];
extern int32_t masVptDimMax_1[1];
extern int32_t masVprDimMax_1[1];
extern int32_t masHelVppDimMax_1[1];
extern int32_t masHelVptDimMax_1[1];
extern int32_t masHelVprDimMax_1[1];

extern int32_t masVtpDimMax_0[1];
extern int32_t masVttDimMax_0[1];
extern int32_t masVtrDimMax_0[1];
extern int32_t masHelVtpDimMax_0[1];
extern int32_t masHelVttDimMax_0[1];
extern int32_t masHelVtrDimMax_0[1];

extern int32_t masVtpDimMax_1[1];
extern int32_t masVttDimMax_1[1];
extern int32_t masVtrDimMax_1[1];
extern int32_t masHelVtpDimMax_1[1];
extern int32_t masHelVttDimMax_1[1];
extern int32_t masHelVtrDimMax_1[1];

extern int32_t masVrpDimMax_0[1];
extern int32_t masVrtDimMax_0[1];
extern int32_t masVrrDimMax_0[1];
extern int32_t masHelVrpDimMax_0[1];
extern int32_t masHelVrtDimMax_0[1];
extern int32_t masHelVrrDimMax_0[1];

extern int32_t masVrpDimMax_1[1];
extern int32_t masVrtDimMax_1[1];
extern int32_t masVrrDimMax_1[1];
extern int32_t masHelVrpDimMax_1[1];
extern int32_t masHelVrtDimMax_1[1];
extern int32_t masHelVrrDimMax_1[1];

extern int32_t masDpDimMax_0[1];
extern int32_t masDtDimMax_0[1];
extern int32_t masDrDimMax_0[1];
extern int32_t masHelDpDimMax_0[1];
extern int32_t masHelDtDimMax_0[1];
extern int32_t masHelDrDimMax_0[1];

extern int32_t masDpDimMax_1[1];
extern int32_t masDtDimMax_1[1];
extern int32_t masDrDimMax_1[1];
extern int32_t masHelDpDimMax_1[1];
extern int32_t masHelDtDimMax_1[1];
extern int32_t masHelDrDimMax_1[1];

extern float * masBppDim_0;
extern float * masBptDim_0;
extern float * masBprDim_0;
extern float * masHelBppDim_0;
extern float * masHelBptDim_0;
extern float * masHelBprDim_0;

extern float * masBtpDim_0;
extern float * masBttDim_0;
extern float * masBtrDim_0;
extern float * masHelBtpDim_0;
extern float * masHelBttDim_0;
extern float * masHelBtrDim_0;

extern float * masBrpDim_0;
extern float * masBrtDim_0;
extern float * masBrrDim_0;
extern float * masHelBrpDim_0;
extern float * masHelBrtDim_0;
extern float * masHelBrrDim_0;

extern float * masVppDim_0;
extern float * masVptDim_0;
extern float * masVprDim_0;
extern float * masHelVppDim_0;
extern float * masHelVptDim_0;
extern float * masHelVprDim_0;

extern float * masVtpDim_0;
extern float * masVttDim_0;
extern float * masVtrDim_0;
extern float * masHelVtpDim_0;
extern float * masHelVttDim_0;
extern float * masHelVtrDim_0;

extern float * masVrpDim_0;
extern float * masVrtDim_0;
extern float * masVrrDim_0;
extern float * masHelVrpDim_0;
extern float * masHelVrtDim_0;
extern float * masHelVrrDim_0;

extern float * masDpDim_0;
extern float * masDtDim_0;
extern float * masDrDim_0;
extern float * masHelDpDim_0;
extern float * masHelDtDim_0;
extern float * masHelDrDim_0;

extern float * masBppDim_1;
extern float * masBptDim_1;
extern float * masBprDim_1;
extern float * masHelBppDim_1;
extern float * masHelBptDim_1;
extern float * masHelBprDim_1;

extern float * masBtpDim_1;
extern float * masBttDim_1;
extern float * masBtrDim_1;
extern float * masHelBtpDim_1;
extern float * masHelBttDim_1;
extern float * masHelBtrDim_1;

extern float * masBrpDim_1;
extern float * masBrtDim_1;
extern float * masBrrDim_1;
extern float * masHelBrpDim_1;
extern float * masHelBrtDim_1;
extern float * masHelBrrDim_1;

extern float * masVppDim_1;
extern float * masVptDim_1;
extern float * masVprDim_1;
extern float * masHelVppDim_1;
extern float * masHelVptDim_1;
extern float * masHelVprDim_1;

extern float * masVtpDim_1;
extern float * masVttDim_1;
extern float * masVtrDim_1;
extern float * masHelVtpDim_1;
extern float * masHelVttDim_1;
extern float * masHelVtrDim_1;

extern float * masVrpDim_1;
extern float * masVrtDim_1;
extern float * masVrrDim_1;
extern float * masHelVrpDim_1;
extern float * masHelVrtDim_1;
extern float * masHelVrrDim_1;

extern float * masDpDim_1;
extern float * masDtDim_1;
extern float * masDrDim_1;
extern float * masHelDpDim_1;
extern float * masHelDtDim_1;
extern float * masHelDrDim_1;

extern float * masBp_0;
extern float * masBt_0;
extern float * masBr_0;
extern float * masVp_0;
extern float * masVt_0;
extern float * masVr_0;
extern float * masD_0;
extern float * masHelBp_0;
extern float * masHelBt_0;
extern float * masHelBr_0;
extern float * masHelVp_0;
extern float * masHelVt_0;
extern float * masHelVr_0;
extern float * masHelD_0;

extern float * masBp_1;
extern float * masBt_1;
extern float * masBr_1;
extern float * masVp_1;
extern float * masVt_1;
extern float * masVr_1;
extern float * masD_1;
extern float * masHelBp_1;
extern float * masHelBt_1;
extern float * masHelBr_1;
extern float * masHelVp_1;
extern float * masHelVt_1;
extern float * masHelVr_1;
extern float * masHelD_1;

/*-- Windows for shared arrays. --*/
extern MPI_Win masBp_0_win;
extern MPI_Win masBt_0_win;
extern MPI_Win masBr_0_win;
extern MPI_Win masVp_0_win;
extern MPI_Win masVt_0_win;
extern MPI_Win masVr_0_win;
extern MPI_Win masD_0_win;
extern MPI_Win masHelBp_0_win;
extern MPI_Win masHelBt_0_win;
extern MPI_Win masHelBr_0_win;
extern MPI_Win masHelVp_0_win;
extern MPI_Win masHelVt_0_win;
extern MPI_Win masHelVr_0_win;
extern MPI_Win masHelD_0_win;

extern MPI_Win masBp_1_win;
extern MPI_Win masBt_1_win;
extern MPI_Win masBr_1_win;
extern MPI_Win masVp_1_win;
extern MPI_Win masVt_1_win;
extern MPI_Win masVr_1_win;
extern MPI_Win masD_1_win;
extern MPI_Win masHelBp_1_win;
extern MPI_Win masHelBt_1_win;
extern MPI_Win masHelBr_1_win;
extern MPI_Win masHelVp_1_win;
extern MPI_Win masHelVt_1_win;
extern MPI_Win masHelVr_1_win;
extern MPI_Win masHelD_1_win;
extern char file_extension[5];

void ERR(int);

void masFetchCouplingInfo(void);

void masFetchFileList(void);

void masGetInterpData( Scalar_t dt );

void masReadData(Index_t fileIndex,
                 float * *masBp, float * *masBt, float * *masBr,
                 float * *masVp, float * *masVt, float * *masVr,
                 float * *masD);

void masHelReadData(Index_t fileIndex,
                    float * *masBp, float * *masBt, float * *masBr,
                    float * *masVp, float * *masVt, float * *masVr,
                    float * *masD);

void masReadFieldIndex(char *state, Index_t fileIndex, Index_t mallocFlag);
void masHelReadFieldIndex(char *state, Index_t fileIndex, Index_t mallocFlag);
void cleanupMPIWindows(char *state);
void cleanupHelMPIWindows(char *state);

void allocateMPIWindows(char *state);
void allocateHelMPIWindows(char *state);

int32_t masReadBppDimTest(Index_t fileIndex);

#endif

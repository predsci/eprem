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
#include <mfhdf.h>
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

extern int32 masDimMin[1];

extern int32 masBppDimMax[1];
extern int32 masBptDimMax[1];
extern int32 masBprDimMax[1];
extern int32 masHelBppDimMax[1];
extern int32 masHelBptDimMax[1];
extern int32 masHelBprDimMax[1];

extern int32 masBtpDimMax[1];
extern int32 masBttDimMax[1];
extern int32 masBtrDimMax[1];
extern int32 masHelBtpDimMax[1];
extern int32 masHelBttDimMax[1];
extern int32 masHelBtrDimMax[1];

extern int32 masBrpDimMax[1];
extern int32 masBrtDimMax[1];
extern int32 masBrrDimMax[1];
extern int32 masHelBrpDimMax[1];
extern int32 masHelBrtDimMax[1];
extern int32 masHelBrrDimMax[1];

extern int32 masVppDimMax[1];
extern int32 masVptDimMax[1];
extern int32 masVprDimMax[1];
extern int32 masHelVppDimMax[1];
extern int32 masHelVptDimMax[1];
extern int32 masHelVprDimMax[1];

extern int32 masVtpDimMax[1];
extern int32 masVttDimMax[1];
extern int32 masVtrDimMax[1];
extern int32 masHelVtpDimMax[1];
extern int32 masHelVttDimMax[1];
extern int32 masHelVtrDimMax[1];

extern int32 masVrpDimMax[1];
extern int32 masVrtDimMax[1];
extern int32 masVrrDimMax[1];
extern int32 masHelVrpDimMax[1];
extern int32 masHelVrtDimMax[1];
extern int32 masHelVrrDimMax[1];

extern int32 masDpDimMax[1];
extern int32 masDtDimMax[1];
extern int32 masDrDimMax[1];
extern int32 masHelDpDimMax[1];
extern int32 masHelDtDimMax[1];
extern int32 masHelDrDimMax[1];

extern float * masBppDim;
extern float * masBptDim;
extern float * masBprDim;
extern float * masHelBppDim;
extern float * masHelBptDim;
extern float * masHelBprDim;

extern float * masBtpDim;
extern float * masBttDim;
extern float * masBtrDim;
extern float * masHelBtpDim;
extern float * masHelBttDim;
extern float * masHelBtrDim;

extern float * masBrpDim;
extern float * masBrtDim;
extern float * masBrrDim;
extern float * masHelBrpDim;
extern float * masHelBrtDim;
extern float * masHelBrrDim;

extern float * masVppDim;
extern float * masVptDim;
extern float * masVprDim;
extern float * masHelVppDim;
extern float * masHelVptDim;
extern float * masHelVprDim;

extern float * masVtpDim;
extern float * masVttDim;
extern float * masVtrDim;
extern float * masHelVtpDim;
extern float * masHelVttDim;
extern float * masHelVtrDim;

extern float * masVrpDim;
extern float * masVrtDim;
extern float * masVrrDim;
extern float * masHelVrpDim;
extern float * masHelVrtDim;
extern float * masHelVrrDim;

extern float * masDpDim;
extern float * masDtDim;
extern float * masDrDim;
extern float * masHelDpDim;
extern float * masHelDtDim;
extern float * masHelDrDim;


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

void ERR(intn);

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

void masReadFieldIndex(void);
void masHelReadFieldIndex(void);

/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/*--*/     void                                                 /*--*/
/*--*/     cleanupMPIWindows(void);                             /*--*/
/*--                                                              --*/
/*--  This function cleans up MPI windows.                        --*/
/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/

#endif

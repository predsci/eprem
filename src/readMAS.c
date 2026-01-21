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

#include <math.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <string.h>
#include <stdlib.h>
#include "readMAS.h"
#include "mpiInit.h"
#include "global.h"
#include "configuration.h"
#include "error.h"
#include "simCore.h"
#include "flow.h"
#include "observerOutput.h"
#include "timers.h"
#include "masIO.h"

Scalar_t *masTime;
Scalar_t *masHelTime;

Scalar_t phiOffset;
Scalar_t phiHelOffset;

Index_t masFileIndex0;
Index_t masFileIndex1;
Index_t masHelFileIndex0;
Index_t masHelFileIndex1;

Index_t masFileIndex_loaded0=-9999;
Index_t masFileIndex_loaded1=-9999;
Index_t masHelFileIndex_loaded0=-9999;
Index_t masHelFileIndex_loaded1=-9999;

Index_t masMallocFlag;
Index_t masEqFileFlag;

int32_t masDimMin[1] = {0};

int32_t masBppDimMax[1];
int32_t masBptDimMax[1];
int32_t masBprDimMax[1];
int32_t masHelBppDimMax[1];
int32_t masHelBptDimMax[1];
int32_t masHelBprDimMax[1];

int32_t masBtpDimMax[1];
int32_t masBttDimMax[1];
int32_t masBtrDimMax[1];
int32_t masHelBtpDimMax[1];
int32_t masHelBttDimMax[1];
int32_t masHelBtrDimMax[1];

int32_t masBrpDimMax[1];
int32_t masBrtDimMax[1];
int32_t masBrrDimMax[1];
int32_t masHelBrpDimMax[1];
int32_t masHelBrtDimMax[1];
int32_t masHelBrrDimMax[1];

int32_t masVppDimMax[1];
int32_t masVptDimMax[1];
int32_t masVprDimMax[1];
int32_t masHelVppDimMax[1];
int32_t masHelVptDimMax[1];
int32_t masHelVprDimMax[1];

int32_t masVtpDimMax[1];
int32_t masVttDimMax[1];
int32_t masVtrDimMax[1];
int32_t masHelVtpDimMax[1];
int32_t masHelVttDimMax[1];
int32_t masHelVtrDimMax[1];

int32_t masVrpDimMax[1];
int32_t masVrtDimMax[1];
int32_t masVrrDimMax[1];
int32_t masHelVrpDimMax[1];
int32_t masHelVrtDimMax[1];
int32_t masHelVrrDimMax[1];

int32_t masDpDimMax[1];
int32_t masDtDimMax[1];
int32_t masDrDimMax[1];
int32_t masHelDpDimMax[1];
int32_t masHelDtDimMax[1];
int32_t masHelDrDimMax[1];

float * masBppDim;
float * masBptDim;
float * masBprDim;
float * masHelBppDim;
float * masHelBptDim;
float * masHelBprDim;

float * masBtpDim;
float * masBttDim;
float * masBtrDim;
float * masHelBtpDim;
float * masHelBttDim;
float * masHelBtrDim;

float * masBrpDim;
float * masBrtDim;
float * masBrrDim;
float * masHelBrpDim;
float * masHelBrtDim;
float * masHelBrrDim;

float * masVppDim;
float * masVptDim;
float * masVprDim;
float * masHelVppDim;
float * masHelVptDim;
float * masHelVprDim;

float * masVtpDim;
float * masVttDim;
float * masVtrDim;
float * masHelVtpDim;
float * masHelVttDim;
float * masHelVtrDim;

float * masVrpDim;
float * masVrtDim;
float * masVrrDim;
float * masHelVrpDim;
float * masHelVrtDim;
float * masHelVrrDim;

float * masDpDim;
float * masDtDim;
float * masDrDim;
float * masHelDpDim;
float * masHelDtDim;
float * masHelDrDim;


float * masBp_0;
float * masBt_0;
float * masBr_0;
float * masVp_0;
float * masVt_0;
float * masVr_0;
float * masD_0;
float * masHelBp_0;
float * masHelBt_0;
float * masHelBr_0;
float * masHelVp_0;
float * masHelVt_0;
float * masHelVr_0;
float * masHelD_0;

float * masBp_1;
float * masBt_1;
float * masBr_1;
float * masVp_1;
float * masVt_1;
float * masVr_1;
float * masD_1;
float * masHelBp_1;
float * masHelBt_1;
float * masHelBr_1;
float * masHelVp_1;
float * masHelVt_1;
float * masHelVr_1;
float * masHelD_1;

MPI_Win masBp_0_win;
MPI_Win masBt_0_win;
MPI_Win masBr_0_win;
MPI_Win masVp_0_win;
MPI_Win masVt_0_win;
MPI_Win masVr_0_win;
MPI_Win masD_0_win;
MPI_Win masHelBp_0_win;
MPI_Win masHelBt_0_win;
MPI_Win masHelBr_0_win;
MPI_Win masHelVp_0_win;
MPI_Win masHelVt_0_win;
MPI_Win masHelVr_0_win;
MPI_Win masHelD_0_win;

MPI_Win masBp_1_win;
MPI_Win masBt_1_win;
MPI_Win masBr_1_win;
MPI_Win masVp_1_win;
MPI_Win masVt_1_win;
MPI_Win masVr_1_win;
MPI_Win masD_1_win;
MPI_Win masHelBp_1_win;
MPI_Win masHelBt_1_win;
MPI_Win masHelBr_1_win;
MPI_Win masHelVp_1_win;
MPI_Win masHelVt_1_win;
MPI_Win masHelVr_1_win;
MPI_Win masHelD_1_win;

int coupleStarted=0;
char file_extension[5];

/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/*--*/     void                                                 /*--*/
/*--*/     ERR(int status)                                     /*--*/
/*--                                                              --*/
/*--  This function checks for an error when loading an SD        --*/
/*--  element.                                                    --*/
/*------------------------------------------------------------------*/
{/*-----------------------------------------------------------------*/
  if (status != 0){
    if (mpi_rank == 0){
      printf("\n\nERROR WITH AN SD! status = %d\n",status);
    }
    exit(1);
  }
}/*-------- END ERR(intn status)  ----------------------------------*/
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/*--*/     void                                                 /*--*/
/*--*/     masFetchCouplingInfo(void)                           /*--*/
/*--                                                              --*/
/*--  This function loads coupling info from the MAS directories  --*/
/*------------------------------------------------------------------*/
{/*-----------------------------------------------------------------*/

  FILE *rfile;

  int max = 10000;
  char line[10000];

  char *name = NULL;
  char *value = NULL;
  char delims[] = ": ";

  char masHelInfoFilename[MAX_STRING_SIZE] = "mas_helio_run_info.txt";
  char masHelInfoFilenameWithPath[MAX_STRING_SIZE];

  // -- coronal coupling --//

  // -- heliospheric coupling --//
  if (config.masHelCouple > 0) {

    // build the path to the info file
    sprintf(masHelInfoFilenameWithPath, "%s%s", config.masHelDirectory, masHelInfoFilename);

    // attempt to open the info file
    rfile = fopen(masHelInfoFilenameWithPath, "r");
    if (rfile==NULL) {
      printf("ERROR - Could not open file \"%s\"\nReverting to defaults\n", masHelInfoFilenameWithPath);
    } else {
      if (mpi_rank == 0) {
        printf("Reading parameters from \"%s\"\n", masHelInfoFilenameWithPath);
      }
    }

    // read each line and check for known parameters
    // this assumes that each line contains `<name>: <value>` or `<name>:<value>`
    while (fgets(line, max, rfile) != NULL) {
      name = strtok(line, delims);
      if (name != NULL) {
        value = strtok(NULL, delims);
        if ((value != NULL) && (strcmp(name, "coronal_helio_phi_offset") == 0)) {
          phiHelOffset = (Scalar_t)atof(value);
        }
      }
    }
    if (mpi_rank == 0) {
      printf("phiHelOffset: %f\n", phiHelOffset);
    }

    // close the file
    fclose(rfile);

    name = NULL;
    value = NULL;

  }

}/*-------- END masFetchCouplingInfo()  ----------------------------*/
/*------------------------------------------------------------------*/

/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/*--*/     void                                                 /*--*/
/*--*/     masFetchFileList(void)                               /*--*/
/*--                                                              --*/
/*--  This function loads the list timesteps and filenames to use --*/
/*--  with MAS files.                                             --*/
/*------------------------------------------------------------------*/
{/*-----------------------------------------------------------------*/

  FILE * rfile;

  int max = 10000;
  char line[10000];

  int nFileLines;

  char *result = NULL;
  char delims[] = " \t\n\r";

  char masTimeFilename[MAX_STRING_SIZE] = "masTime.txt";
  char masTimeFilenameWithPath[MAX_STRING_SIZE];

  Scalar_t initialTime = 0;
  Scalar_t timeVar;

  MPI_Aint N;
  MPI_Aint size;
  int disp_unit;

  // build the path to the time file
  sprintf(masTimeFilenameWithPath, "%s%s", config.masDirectory, masTimeFilename);

  // attempt to open the time file
  rfile = fopen(masTimeFilenameWithPath, "r");
  if (rfile==NULL) {
    printf("ERROR - Could not open file \"%s\"\n", masTimeFilenameWithPath);
    panic("Can't find MAS time list.");
  }

  // read the number of lines in the time file list
  nFileLines = 0;
  while (fgets(line, max, rfile) != NULL) {
    nFileLines++;
  }
  masTime = (Scalar_t *)malloc(sizeof(Scalar_t) * nFileLines);

  // reset the file pointer to the beginning of the file
  rewind(rfile);

  // read in the time file list
  for (int t=0; t<nFileLines; t++) {
    if (fgets(line, max, rfile) != NULL) {
      result = strtok(line, delims);
      if (result != NULL) {
        timeVar = (Scalar_t)atof(result) * MAS_TIME_CONVERT / DAY;
        if (t == 0) {initialTime = timeVar;}
        masTime[t] = config.masStartTime / DAY + (timeVar - initialTime);
      }
      result = strtok(NULL, delims);
    }
  }

  // close the file
  fclose(rfile);

  // store the number of files
  config.masNumFiles = nFileLines;

  result = NULL;

  // if coupling to the heliospheric domain
  if (config.masHelCouple > 0) {

    // build the path to the time file
    sprintf(masTimeFilenameWithPath, "%s%s", config.masHelDirectory, masTimeFilename);

    // attempt to open the time file
    rfile = fopen(masTimeFilenameWithPath, "r");
    if (rfile==NULL) {
      printf("ERROR - Could not open file \"%s\"\n", masTimeFilenameWithPath);
      panic("Can't find MAS Helio time list.");
    }

    // read the number of lines in the time file list
    nFileLines = 0;
    while (fgets(line, max, rfile) != NULL) {
      nFileLines++;
    }
    masHelTime = (Scalar_t *)malloc(sizeof(Scalar_t) * nFileLines);

    // reset the file pointer to the beginning of the file
    rewind(rfile);

    // read in the time file list
    for (int t=0; t<nFileLines; t++) {
      if (fgets(line, max, rfile) != NULL) {
        result = strtok(line, delims);
        if (result != NULL) {
          timeVar = (Scalar_t)atof(result) * MAS_TIME_CONVERT / DAY;
          if (t == 0) {initialTime = timeVar;}
          masHelTime[t] = config.masStartTime / DAY + (timeVar - initialTime);
        }
        result = strtok(NULL, delims);
      }
    }

    // close the file
    fclose(rfile);

    // store the number of files
    config.masHelNumFiles = nFileLines;

    result = NULL;

  }

  // Malloc arrays. (Only malloc once).
  if (masMallocFlag == 0) {

    strncpy(file_extension,".h5",strlen(".h5")+1);

    // The size of the index arrays doesn't change in time in this version.
    masReadFieldIndex();
    // heliospheric coupling
    if (config.masHelCouple > 0)
      masHelReadFieldIndex();

    if(mpi_rank_shared==0){
      N=(int)masBprDimMax[0]*(int)masBptDimMax[0]*(int)masBppDimMax[0];
      MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masBp_0, &masBp_0_win);
      MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masBp_1, &masBp_1_win);
      N=(int)masBtrDimMax[0]*(int)masBttDimMax[0]*(int)masBtpDimMax[0];
      MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masBt_0, &masBt_0_win);
      MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masBt_1, &masBt_1_win);
      N=(int)masBrrDimMax[0]*(int)masBrtDimMax[0]*(int)masBrpDimMax[0];
      MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masBr_0, &masBr_0_win);
      MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masBr_1, &masBr_1_win);
      N=(int)masVprDimMax[0]*(int)masVptDimMax[0]*(int)masVppDimMax[0];
      MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masVp_0, &masVp_0_win);
      MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masVp_1, &masVp_1_win);
      N=(int)masVtrDimMax[0]*(int)masVttDimMax[0]*(int)masVtpDimMax[0];
      MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masVt_0, &masVt_0_win);
      MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masVt_1, &masVt_1_win);
      N=(int)masVrrDimMax[0]*(int)masVrtDimMax[0]*(int)masVrpDimMax[0];
      MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masVr_0, &masVr_0_win);
      MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masVr_1, &masVr_1_win);
      N=(int)masDrDimMax[0]*(int)masDtDimMax[0]*(int)masDpDimMax[0];
      MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masD_0, &masD_0_win);
      MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masD_1, &masD_1_win);

      // heliospheric coupling
      if (config.masHelCouple > 0) {

        N=(int)masHelBprDimMax[0]*(int)masHelBptDimMax[0]*(int)masHelBppDimMax[0];
        MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masHelBp_0, &masHelBp_0_win);
        MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masHelBp_1, &masHelBp_1_win);
        N=(int)masHelBtrDimMax[0]*(int)masHelBttDimMax[0]*(int)masHelBtpDimMax[0];
        MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masHelBt_0, &masHelBt_0_win);
        MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masHelBt_1, &masHelBt_1_win);
        N=(int)masHelBrrDimMax[0]*(int)masHelBrtDimMax[0]*(int)masHelBrpDimMax[0];
        MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masHelBr_0, &masHelBr_0_win);
        MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masHelBr_1, &masHelBr_1_win);
        N=(int)masHelVprDimMax[0]*(int)masHelVptDimMax[0]*(int)masHelVppDimMax[0];
        MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masHelVp_0, &masHelVp_0_win);
        MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masHelVp_1, &masHelVp_1_win);
        N=(int)masHelVtrDimMax[0]*(int)masHelVttDimMax[0]*(int)masHelVtpDimMax[0];
        MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masHelVt_0, &masHelVt_0_win);
        MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masHelVt_1, &masHelVt_1_win);
        N=(int)masHelVrrDimMax[0]*(int)masHelVrtDimMax[0]*(int)masHelVrpDimMax[0];
        MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masHelVr_0, &masHelVr_0_win);
        MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masHelVr_1, &masHelVr_1_win);
        N=(int)masHelDrDimMax[0]*(int)masHelDtDimMax[0]*(int)masHelDpDimMax[0];
        MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masHelD_0, &masHelD_0_win);
        MPI_Win_allocate_shared(N*sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masHelD_1, &masHelD_1_win);

      }

    }else{
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masBp_0, &masBp_0_win);
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masBp_1, &masBp_1_win);
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masBt_0, &masBt_0_win);
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masBt_1, &masBt_1_win);
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masBr_0, &masBr_0_win);
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masBr_1, &masBr_1_win);
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masVp_0, &masVp_0_win);
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masVp_1, &masVp_1_win);
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masVt_0, &masVt_0_win);
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masVt_1, &masVt_1_win);
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masVr_0, &masVr_0_win);
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masVr_1, &masVr_1_win);
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masD_0, &masD_0_win);
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masD_1, &masD_1_win);
      MPI_Win_shared_query(masBp_0_win, 0, &size, &disp_unit, &masBp_0);
      MPI_Win_shared_query(masBt_0_win, 0, &size, &disp_unit, &masBt_0);
      MPI_Win_shared_query(masBr_0_win, 0, &size, &disp_unit, &masBr_0);
      MPI_Win_shared_query(masVp_0_win, 0, &size, &disp_unit, &masVp_0);
      MPI_Win_shared_query(masVt_0_win, 0, &size, &disp_unit, &masVt_0);
      MPI_Win_shared_query(masVr_0_win, 0, &size, &disp_unit, &masVr_0);
      MPI_Win_shared_query(masD_0_win, 0, &size, &disp_unit, &masD_0);
      MPI_Win_shared_query(masBp_1_win, 0, &size, &disp_unit, &masBp_1);
      MPI_Win_shared_query(masBt_1_win, 0, &size, &disp_unit, &masBt_1);
      MPI_Win_shared_query(masBr_1_win, 0, &size, &disp_unit, &masBr_1);
      MPI_Win_shared_query(masVp_1_win, 0, &size, &disp_unit, &masVp_1);
      MPI_Win_shared_query(masVt_1_win, 0, &size, &disp_unit, &masVt_1);
      MPI_Win_shared_query(masVr_1_win, 0, &size, &disp_unit, &masVr_1);
      MPI_Win_shared_query(masD_1_win, 0, &size, &disp_unit, &masD_1);

      // heliospheric coupling
      if (config.masHelCouple > 0) {

        MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masHelBp_0, &masHelBp_0_win);
        MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masHelBp_1, &masHelBp_1_win);
        MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masHelBt_0, &masHelBt_0_win);
        MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masHelBt_1, &masHelBt_1_win);
        MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masHelBr_0, &masHelBr_0_win);
        MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masHelBr_1, &masHelBr_1_win);
        MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masHelVp_0, &masHelVp_0_win);
        MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masHelVp_1, &masHelVp_1_win);
        MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masHelVt_0, &masHelVt_0_win);
        MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masHelVt_1, &masHelVt_1_win);
        MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masHelVr_0, &masHelVr_0_win);
        MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masHelVr_1, &masHelVr_1_win);
        MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masHelD_0, &masHelD_0_win);
        MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masHelD_1, &masHelD_1_win);
        MPI_Win_shared_query(masHelBp_0_win, 0, &size, &disp_unit, &masHelBp_0);
        MPI_Win_shared_query(masHelBt_0_win, 0, &size, &disp_unit, &masHelBt_0);
        MPI_Win_shared_query(masHelBr_0_win, 0, &size, &disp_unit, &masHelBr_0);
        MPI_Win_shared_query(masHelVp_0_win, 0, &size, &disp_unit, &masHelVp_0);
        MPI_Win_shared_query(masHelVt_0_win, 0, &size, &disp_unit, &masHelVt_0);
        MPI_Win_shared_query(masHelVr_0_win, 0, &size, &disp_unit, &masHelVr_0);
        MPI_Win_shared_query(masHelD_0_win, 0, &size, &disp_unit, &masHelD_0);
        MPI_Win_shared_query(masHelBp_1_win, 0, &size, &disp_unit, &masHelBp_1);
        MPI_Win_shared_query(masHelBt_1_win, 0, &size, &disp_unit, &masHelBt_1);
        MPI_Win_shared_query(masHelBr_1_win, 0, &size, &disp_unit, &masHelBr_1);
        MPI_Win_shared_query(masHelVp_1_win, 0, &size, &disp_unit, &masHelVp_1);
        MPI_Win_shared_query(masHelVt_1_win, 0, &size, &disp_unit, &masHelVt_1);
        MPI_Win_shared_query(masHelVr_1_win, 0, &size, &disp_unit, &masHelVr_1);
        MPI_Win_shared_query(masHelD_1_win, 0, &size, &disp_unit, &masHelD_1);

      }

    }
    MPI_Win_lock_all(MPI_MODE_NOCHECK, masBp_0_win);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, masBt_0_win);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, masBr_0_win);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, masVp_0_win);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, masVt_0_win);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, masVr_0_win);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, masD_0_win);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, masBp_1_win);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, masBt_1_win);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, masBr_1_win);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, masVp_1_win);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, masVt_1_win);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, masVr_1_win);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, masD_1_win);

    // heliospheric coupling
    if (config.masHelCouple > 0) {

      MPI_Win_lock_all(MPI_MODE_NOCHECK, masHelBp_0_win);
      MPI_Win_lock_all(MPI_MODE_NOCHECK, masHelBt_0_win);
      MPI_Win_lock_all(MPI_MODE_NOCHECK, masHelBr_0_win);
      MPI_Win_lock_all(MPI_MODE_NOCHECK, masHelVp_0_win);
      MPI_Win_lock_all(MPI_MODE_NOCHECK, masHelVt_0_win);
      MPI_Win_lock_all(MPI_MODE_NOCHECK, masHelVr_0_win);
      MPI_Win_lock_all(MPI_MODE_NOCHECK, masHelD_0_win);
      MPI_Win_lock_all(MPI_MODE_NOCHECK, masHelBp_1_win);
      MPI_Win_lock_all(MPI_MODE_NOCHECK, masHelBt_1_win);
      MPI_Win_lock_all(MPI_MODE_NOCHECK, masHelBr_1_win);
      MPI_Win_lock_all(MPI_MODE_NOCHECK, masHelVp_1_win);
      MPI_Win_lock_all(MPI_MODE_NOCHECK, masHelVt_1_win);
      MPI_Win_lock_all(MPI_MODE_NOCHECK, masHelVr_1_win);
      MPI_Win_lock_all(MPI_MODE_NOCHECK, masHelD_1_win);

    }

    masMallocFlag = 1;

  }


}/*-------- END masFetchFileList()  --------------------------------*/
/*------------------------------------------------------------------*/

/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/*--*/     void                                                 /*--*/
/*--*/     cleanupMPIWindows(void)                              /*--*/
/*--                                                              --*/
/*--  This function cleans up MPI windows.                        --*/
/*------------------------------------------------------------------*/
{/*-----------------------------------------------------------------*/

   if (config.masCouple > 0) {
      MPI_Win_unlock_all(masBp_0_win);
      MPI_Win_unlock_all(masBt_0_win);
      MPI_Win_unlock_all(masBr_0_win);
      MPI_Win_unlock_all(masVp_0_win);
      MPI_Win_unlock_all(masVt_0_win);
      MPI_Win_unlock_all(masVr_0_win);
      MPI_Win_unlock_all(masD_0_win);
      MPI_Win_unlock_all(masBp_1_win);
      MPI_Win_unlock_all(masBt_1_win);
      MPI_Win_unlock_all(masBr_1_win);
      MPI_Win_unlock_all(masVp_1_win);
      MPI_Win_unlock_all(masVt_1_win);
      MPI_Win_unlock_all(masVr_1_win);
      MPI_Win_unlock_all(masD_1_win);

      MPI_Win_free(&masBp_0_win);
      MPI_Win_free(&masBt_0_win);
      MPI_Win_free(&masBr_0_win);
      MPI_Win_free(&masVp_0_win);
      MPI_Win_free(&masVt_0_win);
      MPI_Win_free(&masVr_0_win);
      MPI_Win_free(&masD_0_win);
      MPI_Win_free(&masBp_1_win);
      MPI_Win_free(&masBt_1_win);
      MPI_Win_free(&masBr_1_win);
      MPI_Win_free(&masVp_1_win);
      MPI_Win_free(&masVt_1_win);
      MPI_Win_free(&masVr_1_win);
      MPI_Win_free(&masD_1_win);
  }

  // heliospheric coupling
  if (config.masHelCouple > 0) {

    MPI_Win_unlock_all(masHelBp_0_win);
    MPI_Win_unlock_all(masHelBt_0_win);
    MPI_Win_unlock_all(masHelBr_0_win);
    MPI_Win_unlock_all(masHelVp_0_win);
    MPI_Win_unlock_all(masHelVt_0_win);
    MPI_Win_unlock_all(masHelVr_0_win);
    MPI_Win_unlock_all(masHelD_0_win);
    MPI_Win_unlock_all(masHelBp_1_win);
    MPI_Win_unlock_all(masHelBt_1_win);
    MPI_Win_unlock_all(masHelBr_1_win);
    MPI_Win_unlock_all(masHelVp_1_win);
    MPI_Win_unlock_all(masHelVt_1_win);
    MPI_Win_unlock_all(masHelVr_1_win);
    MPI_Win_unlock_all(masHelD_1_win);

    MPI_Win_free(&masHelBp_0_win);
    MPI_Win_free(&masHelBt_0_win);
    MPI_Win_free(&masHelBr_0_win);
    MPI_Win_free(&masHelVp_0_win);
    MPI_Win_free(&masHelVt_0_win);
    MPI_Win_free(&masHelVr_0_win);
    MPI_Win_free(&masHelD_0_win);
    MPI_Win_free(&masHelBp_1_win);
    MPI_Win_free(&masHelBt_1_win);
    MPI_Win_free(&masHelBr_1_win);
    MPI_Win_free(&masHelVp_1_win);
    MPI_Win_free(&masHelVt_1_win);
    MPI_Win_free(&masHelVr_1_win);
    MPI_Win_free(&masHelD_1_win);

  }

}/*-------- END cleanupMPIWindows()  -------------------------------*/
/*------------------------------------------------------------------*/

/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/*--*/     void                                                 /*--*/
/*--*/     masGetInterpData ( Scalar_t dt )                     /*--*/
/*--                                                              --*/
/*-- This function checks if the simulation time is right to read --*/
/*-- or copy MAS data, and does so if needed.                     --*/
/*-- It also sets the interpolation factors 's_cor' and 's_hel'   --*/
/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
{

  int i,N;
  int need_sync_0,need_sync_1;
  Scalar_t time_interp;

//
// Set the time we are interpolating to (the current time plus the time step).
//

  time_interp = t_global + dt;

//
// Check if MAS is being used yet.
//
  if ((t_global < masTime[0]) && (config.masSteadyState == 0)){
    mhdGridStatus = MHD_DEFAULT;
    return;
  }
  else{
    mhdGridStatus = MHD_MAS;
  }

//
// Find the two MAS files that bound the current time_interp, and set
// interpolation factors.
//

  if (time_interp <= masTime[0]){

    masFileIndex0 = 0;
    masFileIndex1 = 0;

    s_cor = 0.0;

  } else if (time_interp >= masTime[config.masNumFiles-1]){

    masFileIndex0 = config.masNumFiles-1;
    masFileIndex1 = config.masNumFiles-1;

    s_cor = 0.0;

  } else {

    for (i=1; i<config.masNumFiles; i++){
      if (masTime[i] > time_interp ){
         masFileIndex0 = i-1;
         masFileIndex1 = i;
         break;
      }
    }

    s_cor = ( time_interp - masTime[masFileIndex0] ) /
            ( masTime[masFileIndex1] - masTime[masFileIndex0]);

  }

//
// Read or copy MAS data as needed.
//

  if (masFileIndex0 != masFileIndex_loaded0){ // if state0 needs to be updated
    MPI_Barrier(comm_shared);
    if (masFileIndex0 == masFileIndex_loaded1){ // if state0 used is already stored previously in state1
      if(mpi_rank_shared==0) {
        N=(int)masBprDimMax[0]*(int)masBptDimMax[0]*(int)masBppDimMax[0];
        memcpy(&masBp_0[0],&masBp_1[0],N*sizeof(float));
        N=(int)masBtrDimMax[0]*(int)masBttDimMax[0]*(int)masBtpDimMax[0];
        memcpy(&masBt_0[0],&masBt_1[0],N*sizeof(float));
        N=(int)masBrrDimMax[0]*(int)masBrtDimMax[0]*(int)masBrpDimMax[0];
        memcpy(&masBr_0[0],&masBr_1[0],N*sizeof(float));
        N=(int)masVprDimMax[0]*(int)masVptDimMax[0]*(int)masVppDimMax[0];
        memcpy(&masVp_0[0],&masVp_1[0],N*sizeof(float));
        N=(int)masVtrDimMax[0]*(int)masVttDimMax[0]*(int)masVtpDimMax[0];
        memcpy(&masVt_0[0],&masVt_1[0],N*sizeof(float));
        N=(int)masVrrDimMax[0]*(int)masVrtDimMax[0]*(int)masVrpDimMax[0];
        memcpy(&masVr_0[0],&masVr_1[0],N*sizeof(float));
        N=(int)masDrDimMax[0]*(int)masDtDimMax[0]*(int)masDpDimMax[0];
        memcpy(&masD_0[0],&masD_1[0],N*sizeof(float));
      }
    } else {
      if(mpi_rank_shared==0) {
        masReadData(masFileIndex0,
                &masBp_0, &masBt_0, &masBr_0,
                &masVp_0, &masVt_0, &masVr_0,
                &masD_0);
      }
    }
    masFileIndex_loaded0 = masFileIndex0;
    need_sync_0 = 1;
  } else {
    need_sync_0 = 0;
  }

  if (masFileIndex1 != masFileIndex_loaded1){
    MPI_Barrier(comm_shared);
    if(mpi_rank_shared==0) {
      masReadData(masFileIndex1,
                &masBp_1, &masBt_1, &masBr_1,
                &masVp_1, &masVt_1, &masVr_1,
                &masD_1);
    }
    masFileIndex_loaded1 = masFileIndex1;
    need_sync_1 = 1;
  } else {
    need_sync_1 = 0;
  }

  if (need_sync_0 == 1){
    MPI_Barrier(comm_shared);
    MPI_Win_sync(masBp_0_win);
    MPI_Win_sync(masBt_0_win);
    MPI_Win_sync(masBr_0_win);
    MPI_Win_sync(masVp_0_win);
    MPI_Win_sync(masVt_0_win);
    MPI_Win_sync(masVr_0_win);
    MPI_Win_sync(masD_0_win);
  }

  if (need_sync_1 == 1){
    MPI_Barrier(comm_shared);
    MPI_Win_sync(masBp_1_win);
    MPI_Win_sync(masBt_1_win);
    MPI_Win_sync(masBr_1_win);
    MPI_Win_sync(masVp_1_win);
    MPI_Win_sync(masVt_1_win);
    MPI_Win_sync(masVr_1_win);
    MPI_Win_sync(masD_1_win);
  }

  if (need_sync_0 + need_sync_1 > 0) MPI_Barrier(comm_shared);
//
// Now do everything with the heliosphere.
//
  if (config.masHelCouple > 0) {

//
// Find the two MAS HELIO files that bound the current time_interp, and set
// interpolation factors.
//
    if (time_interp < masHelTime[0]){

      masHelFileIndex0 = 0;
      masHelFileIndex1 = 0;

      s_hel = 0.0;

    } else if (time_interp > masHelTime[config.masHelNumFiles-1]){

      masHelFileIndex0 = config.masHelNumFiles-1;
      masHelFileIndex1 = config.masHelNumFiles-1;

      s_hel = 0.0;

    } else {

      for (i=1; i<config.masHelNumFiles; i++){
        if (masHelTime[i] > time_interp){
          masHelFileIndex0 = i-1;
          masHelFileIndex1 = i;
          break;
        }
      }

      s_hel = ( time_interp - masHelTime[masHelFileIndex0] ) /
              ( masHelTime[masHelFileIndex1] - masHelTime[masHelFileIndex0]);

    }

//
// Read or copy MAS HELIO data as needed.
//
    if (masHelFileIndex0 != masHelFileIndex_loaded0){// If new data needs to be read into state 0
      MPI_Barrier(comm_shared);
      if (masHelFileIndex0 == masHelFileIndex_loaded1)// If data already exists in previous state 1
      {
        if(mpi_rank_shared==0) {
          N=(int)masHelBprDimMax[0]*(int)masHelBptDimMax[0]*(int)masHelBppDimMax[0];
          memcpy(&masHelBp_0[0],&masHelBp_1[0],N*sizeof(float));
          N=(int)masHelBtrDimMax[0]*(int)masHelBttDimMax[0]*(int)masHelBtpDimMax[0];
          memcpy(&masHelBt_0[0],&masHelBt_1[0],N*sizeof(float));
          N=(int)masHelBrrDimMax[0]*(int)masHelBrtDimMax[0]*(int)masHelBrpDimMax[0];
          memcpy(&masHelBr_0[0],&masHelBr_1[0],N*sizeof(float));
          N=(int)masHelVprDimMax[0]*(int)masHelVptDimMax[0]*(int)masHelVppDimMax[0];
          memcpy(&masHelVp_0[0],&masHelVp_1[0],N*sizeof(float));
          N=(int)masHelVtrDimMax[0]*(int)masHelVttDimMax[0]*(int)masHelVtpDimMax[0];
          memcpy(&masHelVt_0[0],&masHelVt_1[0],N*sizeof(float));
          N=(int)masHelVrrDimMax[0]*(int)masHelVrtDimMax[0]*(int)masHelVrpDimMax[0];
          memcpy(&masHelVr_0[0],&masHelVr_1[0],N*sizeof(float));
          N=(int)masHelDrDimMax[0]*(int)masHelDtDimMax[0]*(int)masHelDpDimMax[0];
          memcpy(&masHelD_0[0],&masHelD_1[0],N*sizeof(float));
        }
      }
      else {// Read in new data
        if(mpi_rank_shared==0) {
          masHelReadData(masHelFileIndex0,
                  &masHelBp_0, &masHelBt_0, &masHelBr_0,
                  &masHelVp_0, &masHelVt_0, &masHelVr_0,
                  &masHelD_0);
        }
      }
      masHelFileIndex_loaded0 = masHelFileIndex0;
      need_sync_0 = 1;
    } else {// Nothing changed
      need_sync_0 = 0;
    }

    if (masHelFileIndex1 != masHelFileIndex_loaded1){ //  Check if state1 needs to be read in
      MPI_Barrier(comm_shared);
      if(mpi_rank_shared==0) {
        masHelReadData(masHelFileIndex1,
                  &masHelBp_1, &masHelBt_1, &masHelBr_1,
                  &masHelVp_1, &masHelVt_1, &masHelVr_1,
                  &masHelD_1);
      }
      masHelFileIndex_loaded1 = masHelFileIndex1;
      need_sync_1 = 1;
    } else {
      need_sync_1 = 0;
    }

    if (need_sync_0 == 1){
      MPI_Barrier(comm_shared);
      MPI_Win_sync(masHelBp_0_win);
      MPI_Win_sync(masHelBt_0_win);
      MPI_Win_sync(masHelBr_0_win);
      MPI_Win_sync(masHelVp_0_win);
      MPI_Win_sync(masHelVt_0_win);
      MPI_Win_sync(masHelVr_0_win);
      MPI_Win_sync(masHelD_0_win);
    }

    if (need_sync_1 == 1){
      MPI_Barrier(comm_shared);
      MPI_Win_sync(masHelBp_1_win);
      MPI_Win_sync(masHelBt_1_win);
      MPI_Win_sync(masHelBr_1_win);
      MPI_Win_sync(masHelVp_1_win);
      MPI_Win_sync(masHelVt_1_win);
      MPI_Win_sync(masHelVr_1_win);
      MPI_Win_sync(masHelD_1_win);
    }

    if (need_sync_0 + need_sync_1 > 0) MPI_Barrier(comm_shared);

  }
//
// Tell everyone what we just did.
//
  if (mpi_rank==0) {
    printf("  --> MAS(C) Couple: Time: %14.8e  Time+Dt: %14.8e  s_cor: %10.8f\n",
          t_global,
          time_interp,s_cor);
    printf("  -->                idx0: %03d  idx1: %03d  masTime0: %14.8e  masTime1: %14.8e\n",
          masFileIndex0+1,
          masFileIndex1+1,
          masTime[masFileIndex0],
          masTime[masFileIndex1]);
    if (config.masHelCouple > 0){
      printf("  --> MAS(H) Couple: Time: %14.8e  Time+Dt: %14.8e  s_hel: %10.8f\n",
          t_global,
          time_interp,s_hel);
      printf("  -->                idx0: %03d  idx1: %03d  masTime0: %14.8e  masTime1: %14.8e\n",
          masHelFileIndex0+1,
          masHelFileIndex1+1,
          masHelTime[masHelFileIndex0],
          masHelTime[masHelFileIndex1]);
    }
  }

}
/*----------------- END masGetInterpData(dt)  ------------------------*/
/*--------------------------------------------------------------------*/

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/
/*--*/     void                                                   /*--*/
/*--*/     masReadData(Index_t fileIndex,                         /*--*/
/*--*/            float *masBp[], float *masBt[], float *masBr[], /*--*/
/*--*/            float *masVp[], float *masVt[], float *masVr[], /*--*/
/*--*/            float *masD[])                                  /*--*/
/*--*/                                                            /*--*/
/*--                                                                --*/
/*--This function reads the MAS data from a HDF file.              --*/
/*--------------------------------------------------------------------*/
{/*-------------------------------------------------------------------*/
  char fileNames[7][MAX_STRING_SIZE];

  double timer_tmp = 0;

  timer_tmp = MPI_Wtime();

  if (config.masDigits == 3) {

    sprintf(fileNames[0], "%sbp%03d%s", config.masDirectory, fileIndex + 1, file_extension);
    sprintf(fileNames[1], "%sbt%03d%s", config.masDirectory, fileIndex + 1, file_extension);
    sprintf(fileNames[2], "%sbr%03d%s", config.masDirectory, fileIndex + 1, file_extension);
    sprintf(fileNames[3], "%svp%03d%s", config.masDirectory, fileIndex + 1, file_extension);
    sprintf(fileNames[4], "%svt%03d%s", config.masDirectory, fileIndex + 1, file_extension);
    sprintf(fileNames[5], "%svr%03d%s", config.masDirectory, fileIndex + 1, file_extension);
    sprintf(fileNames[6], "%srho%03d%s", config.masDirectory, fileIndex + 1, file_extension);

  } else {

    sprintf(fileNames[0], "%sbp%06d%s", config.masDirectory, fileIndex + 1, file_extension);
    sprintf(fileNames[1], "%sbt%06d%s", config.masDirectory, fileIndex + 1, file_extension);
    sprintf(fileNames[2], "%sbr%06d%s", config.masDirectory, fileIndex + 1, file_extension);
    sprintf(fileNames[3], "%svp%06d%s", config.masDirectory, fileIndex + 1, file_extension);
    sprintf(fileNames[4], "%svt%06d%s", config.masDirectory, fileIndex + 1, file_extension);
    sprintf(fileNames[5], "%svr%06d%s", config.masDirectory, fileIndex + 1, file_extension);
    sprintf(fileNames[6], "%srho%06d%s", config.masDirectory, fileIndex + 1, file_extension);

  }

  // reading in Bp
  masReadDatafromFile(fileNames[0], masBp);

  // reading in Bt
  masReadDatafromFile(fileNames[1], masBt);

  // reading in Br
  masReadDatafromFile(fileNames[2], masBr);

  // reading in Vp
  masReadDatafromFile(fileNames[3], masVp);

  // reading in Vt
  masReadDatafromFile(fileNames[4], masVt);

  // reading in Vr
  masReadDatafromFile(fileNames[5], masVr);

  // reading in D
  masReadDatafromFile(fileNames[6], masD);

  if (mpi_rank == 0) printf("  --> IO MAS: Coronal sequence %03d read.\n",fileIndex+1);

  timer_mas_io = timer_mas_io + (MPI_Wtime() - timer_tmp);

}/*-------- END masReadData()  -----------------------*/
/*------------------------------------------------------------------*/


/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/
/*--*/     void                                                   /*--*/
/*--*/     masHelReadData(Index_t fileIndex,                      /*--*/
/*--*/            float *masBp[], float *masBt[], float *masBr[], /*--*/
/*--*/            float *masVp[], float *masVt[], float *masVr[], /*--*/
/*--*/            float *masD[])                                  /*--*/
/*--*/                                                            /*--*/
/*--                                                                --*/
/*--This function checks if the simulation time is right to load    --*/
/*--another timestep file from the MAS HEL data, and does if so.    --*/
/*--------------------------------------------------------------------*/
{/*-------------------------------------------------------------------*/

  char fileNames[7][MAX_STRING_SIZE];

  double timer_tmp = 0;

  timer_tmp = MPI_Wtime();

  if (config.masHelDigits == 3) {

    sprintf(fileNames[0], "%sbp%03d%s", config.masHelDirectory, fileIndex + 1, file_extension);
    sprintf(fileNames[1], "%sbt%03d%s", config.masHelDirectory, fileIndex + 1, file_extension);
    sprintf(fileNames[2], "%sbr%03d%s", config.masHelDirectory, fileIndex + 1, file_extension);
    sprintf(fileNames[3], "%svp%03d%s", config.masHelDirectory, fileIndex + 1, file_extension);
    sprintf(fileNames[4], "%svt%03d%s", config.masHelDirectory, fileIndex + 1, file_extension);
    sprintf(fileNames[5], "%svr%03d%s", config.masHelDirectory, fileIndex + 1, file_extension);
    sprintf(fileNames[6], "%srho%03d%s", config.masHelDirectory, fileIndex + 1, file_extension);
  }
  else
  {
    sprintf(fileNames[0], "%sbp%06d%s", config.masHelDirectory, fileIndex + 1, file_extension);
    sprintf(fileNames[1], "%sbt%06d%s", config.masHelDirectory, fileIndex + 1, file_extension);
    sprintf(fileNames[2], "%sbr%06d%s", config.masHelDirectory, fileIndex + 1, file_extension);
    sprintf(fileNames[3], "%svp%06d%s", config.masHelDirectory, fileIndex + 1, file_extension);
    sprintf(fileNames[4], "%svt%06d%s", config.masHelDirectory, fileIndex + 1, file_extension);
    sprintf(fileNames[5], "%svr%06d%s", config.masHelDirectory, fileIndex + 1, file_extension);
    sprintf(fileNames[6], "%srho%06d%s", config.masHelDirectory, fileIndex + 1, file_extension);
  }

  // reading in Bp
  masReadDatafromFile(fileNames[0], masBp);

  // reading in Bt
  masReadDatafromFile(fileNames[1], masBt);

  // reading in Br
  masReadDatafromFile(fileNames[2], masBr);

  // reading in Vp
  masReadDatafromFile(fileNames[3], masVp);

  // reading in Vt
  masReadDatafromFile(fileNames[4], masVt);

  // reading in Vr
  masReadDatafromFile(fileNames[5], masVr);

  // reading in D
  masReadDatafromFile(fileNames[6], masD);

  if (mpi_rank == 0) printf("  --> IO MAS: Helio sequence %03d read.\n",fileIndex+1);

  timer_mas_io = timer_mas_io + (MPI_Wtime() - timer_tmp);

}/*-------- END masHelReadData()  ----------------------------------*/
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/*--*/     void                                                 /*--*/
/*--*/     masReadFieldIndex(void)                              /*--*/
/*--                                                              --*/
/*------------------------------------------------------------------*/
{/*-----------------------------------------------------------------*/

  char fileNames[7][MAX_STRING_SIZE];

  Scalar_t rMin, rMax, rTemp;

  double timer_tmp=0;

  timer_tmp = MPI_Wtime();

  if (config.masDigits == 3) {

    sprintf(fileNames[0], "%sbp%03d%s", config.masDirectory, 1, file_extension);
    sprintf(fileNames[1], "%sbt%03d%s", config.masDirectory, 1, file_extension);
    sprintf(fileNames[2], "%sbr%03d%s", config.masDirectory, 1, file_extension);
    sprintf(fileNames[3], "%svp%03d%s", config.masDirectory, 1, file_extension);
    sprintf(fileNames[4], "%svt%03d%s", config.masDirectory, 1, file_extension);
    sprintf(fileNames[5], "%svr%03d%s", config.masDirectory, 1, file_extension);
    sprintf(fileNames[6], "%srho%03d%s", config.masDirectory, 1, file_extension);

  } else {

    sprintf(fileNames[0], "%sbp%06d%s", config.masDirectory, 1, file_extension);
    sprintf(fileNames[1], "%sbt%06d%s", config.masDirectory, 1, file_extension);
    sprintf(fileNames[2], "%sbr%06d%s", config.masDirectory, 1, file_extension);
    sprintf(fileNames[3], "%svp%06d%s", config.masDirectory, 1, file_extension);
    sprintf(fileNames[4], "%svt%06d%s", config.masDirectory, 1, file_extension);
    sprintf(fileNames[5], "%svr%06d%s", config.masDirectory, 1, file_extension);
    sprintf(fileNames[6], "%srho%06d%s", config.masDirectory, 1, file_extension);

  }

  // Read in dimensions --> Allocate the mesh storage --> Read the mesh

    // masBppDim
    masReadMeshDimensions(fileNames[0], "dim3", &masBppDimMax[0]);
    masBppDim = (float *)malloc(sizeof(float) * (int)(masBppDimMax[0])); 
    masReadMesh(fileNames[0], "dim3", &masBppDim);

    // masBptDim
    masReadMeshDimensions(fileNames[0], "dim2", &masBptDimMax[0]);
    masBptDim = (float *)malloc(sizeof(float) * (int)(masBptDimMax[0]));
    masReadMesh(fileNames[0], "dim2", &masBptDim);

    // masBprDim
    masReadMeshDimensions(fileNames[0], "dim1", &masBprDimMax[0]);
    masBprDim = (float *)malloc(sizeof(float) * (int)(masBprDimMax[0]));
    masReadMesh(fileNames[0], "dim1", &masBprDim);

    // grab the min and max for the r
    rMin = masBprDim[0];
    rMax = masBprDim[masBprDimMax[0] - 1];

    // masBtpDim
    masReadMeshDimensions(fileNames[1], "dim3", &masBtpDimMax[0]);
    masBtpDim = (float *)malloc(sizeof(float) * (int)(masBtpDimMax[0]));
    masReadMesh(fileNames[1], "dim3", &masBtpDim);

    // masBttDim
    masReadMeshDimensions(fileNames[1], "dim2", &masBttDimMax[0]);
    masBttDim = (float *)malloc(sizeof(float) * (int)(masBttDimMax[0]));
    masReadMesh(fileNames[1], "dim2", &masBttDim);

    // masBtrDim
    masReadMeshDimensions(fileNames[1], "dim1", &masBtrDimMax[0]);
    masBtrDim = (float *)malloc(sizeof(float) * (int)(masBtrDimMax[0]));
    masReadMesh(fileNames[1], "dim1", &masBtrDim);

    // grab the min and max for the r
    rTemp = masBtrDim[0];
    if (rTemp > rMin)
      rMin = rTemp;

    rTemp = masBtrDim[masBtrDimMax[0] - 1];
    if (rTemp < rMax)
      rMax = rTemp;

    // masBrpDim
    masReadMeshDimensions(fileNames[2], "dim3", &masBrpDimMax[0]);
    masBrpDim = (float *)malloc(sizeof(float) * (int)(masBrpDimMax[0]));
    masReadMesh(fileNames[2], "dim3", &masBrpDim);

    // masBrtDim
    masReadMeshDimensions(fileNames[2], "dim2", &masBrtDimMax[0]);
    masBrtDim = (float *)malloc(sizeof(float) * (int)(masBrtDimMax[0]));
    masReadMesh(fileNames[2], "dim2", &masBrtDim);

    // masBrrDim
    masReadMeshDimensions(fileNames[2], "dim1", &masBrrDimMax[0]);
    masBrrDim = (float *)malloc(sizeof(float) * (int)(masBrrDimMax[0]));
    masReadMesh(fileNames[2], "dim1", &masBrrDim);

    // grab the min and max for the r
    rTemp = masBrrDim[0];
    if (rTemp > rMin)
      rMin = rTemp;

    rTemp = masBrrDim[masBrrDimMax[0] - 1];
    if (rTemp < rMax)
      rMax = rTemp;

    // masVppDim
    masReadMeshDimensions(fileNames[3], "dim3", &masVppDimMax[0]);
    masVppDim = (float *)malloc(sizeof(float) * (int)(masVppDimMax[0]));
    masReadMesh(fileNames[3], "dim3", &masVppDim);

    // masVptDim
    masReadMeshDimensions(fileNames[3], "dim2", &masVptDimMax[0]);
    masVptDim = (float *)malloc(sizeof(float) * (int)(masVptDimMax[0]));
    masReadMesh(fileNames[3], "dim2", &masVptDim);

    // masVprDim
    masReadMeshDimensions(fileNames[3], "dim1", &masVprDimMax[0]);
    masVprDim = (float *)malloc(sizeof(float) * (int)(masVprDimMax[0]));
    masReadMesh(fileNames[3], "dim1", &masVprDim);

    // grab the min and max for the r
    rTemp = masVprDim[0];
    if (rTemp > rMin)
      rMin = rTemp;

    rTemp = masVprDim[masVprDimMax[0] - 1];
    if (rTemp < rMax)
      rMax = rTemp;

    // masVtpDim
    masReadMeshDimensions(fileNames[4], "dim3", &masVtpDimMax[0]);
    masVtpDim = (float *)malloc(sizeof(float) * (int)(masVtpDimMax[0]));
    masReadMesh(fileNames[4], "dim3", &masVtpDim);

    // masVttDim
    masReadMeshDimensions(fileNames[4], "dim2", &masVttDimMax[0]);
    masVttDim = (float *)malloc(sizeof(float) * (int)(masVttDimMax[0]));
    masReadMesh(fileNames[4], "dim2", &masVttDim);

    // masVtrDim
    masReadMeshDimensions(fileNames[4], "dim1", &masVtrDimMax[0]);
    masVtrDim = (float *)malloc(sizeof(float) * (int)(masVtrDimMax[0]));
    masReadMesh(fileNames[4], "dim1", &masVtrDim);

    // grab the min and max for the r
    rTemp = masVtrDim[0];
    if (rTemp > rMin)
      rMin = rTemp;

    rTemp = masVtrDim[masVtrDimMax[0] - 1];
    if (rTemp < rMax)
      rMax = rTemp;

    // masVrpDim
    masReadMeshDimensions(fileNames[5], "dim3", &masVrpDimMax[0]);
    masVrpDim = (float *)malloc(sizeof(float) * (int)(masVrpDimMax[0]));
    masReadMesh(fileNames[5], "dim3", &masVrpDim);

    // masVrtDim
    masReadMeshDimensions(fileNames[5], "dim2", &masVrtDimMax[0]);
    masVrtDim = (float *)malloc(sizeof(float) * (int)(masVrtDimMax[0]));
    masReadMesh(fileNames[5], "dim2", &masVrtDim);

    // masVrrDim
    masReadMeshDimensions(fileNames[5], "dim1", &masVrrDimMax[0]);
    masVrrDim = (float *)malloc(sizeof(float) * (int)(masVrrDimMax[0]));
    masReadMesh(fileNames[5], "dim1", &masVrrDim);

    // grab the min and max for the r
    rTemp = masVrrDim[0];
    if (rTemp > rMin)
      rMin = rTemp;

    rTemp = masVrrDim[masVrrDimMax[0] - 1];
    if (rTemp < rMax)
      rMax = rTemp;

    // masDpDim
    masReadMeshDimensions(fileNames[6], "dim3", &masDpDimMax[0]);
    masDpDim = (float *)malloc(sizeof(float) * (int)(masDpDimMax[0]));
    masReadMesh(fileNames[6], "dim3", &masDpDim);

    // masDtDim
    masReadMeshDimensions(fileNames[6], "dim2", &masDtDimMax[0]);
    masDtDim = (float *)malloc(sizeof(float) * (int)(masDtDimMax[0]));
    masReadMesh(fileNames[6], "dim2", &masDtDim);

    // masDrDim
    masReadMeshDimensions(fileNames[6], "dim1", &masDrDimMax[0]);
    masDrDim = (float *)malloc(sizeof(float) * (int)(masDrDimMax[0]));
    masReadMesh(fileNames[6], "dim1", &masDrDim);

    // grab the min and max for the r
    rTemp = masDrDim[0];
    if (rTemp > rMin)
      rMin = rTemp;

    rTemp = masDrDim[masDrDimMax[0] - 1];
    if (rTemp < rMax)
      rMax = rTemp;

  // set rScale, and masRadialMin/Max
  config.rScale = rMin * RSAU;
  config.masRadialMin = rMin;
  config.masRadialMax = rMax;

  timer_mas_io = timer_mas_io + (MPI_Wtime() - timer_tmp);

}/*-------- END masReadFieldIndex()  -----------------------*/
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/*--*/     void                                                 /*--*/
/*--*/     masHelReadFieldIndex(void)                           /*--*/
/*--                                                              --*/
/*------------------------------------------------------------------*/
{/*-----------------------------------------------------------------*/

  char fileNames[7][MAX_STRING_SIZE];
  Scalar_t rMin, rMax, rTemp;

  double timer_tmp = 0;

  timer_tmp = MPI_Wtime();

  if (config.masHelDigits == 3)
  {
    sprintf(fileNames[0], "%sbp%03d%s", config.masHelDirectory, 1, file_extension);
    sprintf(fileNames[1], "%sbt%03d%s", config.masHelDirectory, 1, file_extension);
    sprintf(fileNames[2], "%sbr%03d%s", config.masHelDirectory, 1, file_extension);
    sprintf(fileNames[3], "%svp%03d%s", config.masHelDirectory, 1, file_extension);
    sprintf(fileNames[4], "%svt%03d%s", config.masHelDirectory, 1, file_extension);
    sprintf(fileNames[5], "%svr%03d%s", config.masHelDirectory, 1, file_extension);
    sprintf(fileNames[6], "%srho%03d%s", config.masHelDirectory,1, file_extension);
  }
  else
  {
    sprintf(fileNames[0], "%sbp%06d%s", config.masHelDirectory,  1, file_extension);
    sprintf(fileNames[1], "%sbt%06d%s", config.masHelDirectory,  1, file_extension);
    sprintf(fileNames[2], "%sbr%06d%s", config.masHelDirectory,  1, file_extension);
    sprintf(fileNames[3], "%svp%06d%s", config.masHelDirectory,  1, file_extension);
    sprintf(fileNames[4], "%svt%06d%s", config.masHelDirectory,  1, file_extension);
    sprintf(fileNames[5], "%svr%06d%s", config.masHelDirectory,  1, file_extension);
    sprintf(fileNames[6], "%srho%06d%s", config.masHelDirectory, 1, file_extension);
  }

    // masHelBpp
    masReadMeshDimensions(fileNames[0], "dim3", &masHelBppDimMax[0]);
    masHelBppDim = (float *)malloc(sizeof(float) * (int)(masHelBppDimMax[0]));
    masReadMesh(fileNames[0], "dim3", &masHelBppDim);

    // masHelBptDim
    masReadMeshDimensions(fileNames[0], "dim2", &masHelBptDimMax[0]);
    masHelBptDim = (float *)malloc(sizeof(float) * (int)(masHelBptDimMax[0]));
    masReadMesh(fileNames[0], "dim2", &masHelBptDim);

    // masHelBprDim
    masReadMeshDimensions(fileNames[0], "dim1", &masHelBprDimMax[0]);
    masHelBprDim = (float *)malloc(sizeof(float) * (int)(masHelBprDimMax[0]));
    masReadMesh(fileNames[0], "dim1", &masHelBprDim);

    // grab the min and max for the r
    rMin = masHelBprDim[0];
    rMax = masHelBprDim[masHelBprDimMax[0] - 1];

    // masHelBtpDim
    masReadMeshDimensions(fileNames[1], "dim3", &masHelBtpDimMax[0]);
    masHelBtpDim = (float *)malloc(sizeof(float) * (int)(masHelBtpDimMax[0]));
    masReadMesh(fileNames[1], "dim3", &masHelBtpDim);

    // masHelBttDim
    masReadMeshDimensions(fileNames[1], "dim2", &masHelBttDimMax[0]);
    masHelBttDim = (float *)malloc(sizeof(float) * (int)(masHelBttDimMax[0]));
    masReadMesh(fileNames[1], "dim2", &masHelBttDim);

    // masHelBtrDim
    masReadMeshDimensions(fileNames[1], "dim1", &masHelBtrDimMax[0]);
    masHelBtrDim = (float *)malloc(sizeof(float) * (int)(masHelBtrDimMax[0]));
    masReadMesh(fileNames[1], "dim1", &masHelBtrDim);

    // grab the min and max for the r
    rTemp = masHelBtrDim[0];
    if (rTemp > rMin)
      rMin = rTemp;

    rTemp = masHelBtrDim[masHelBtrDimMax[0] - 1];
    if (rTemp < rMax)
      rMax = rTemp;

    // masHelBrpDim
    masReadMeshDimensions(fileNames[2], "dim3", &masHelBrpDimMax[0]);
    masHelBrpDim = (float *)malloc(sizeof(float) * (int)(masHelBrpDimMax[0]));
    masReadMesh(fileNames[2], "dim3", &masHelBrpDim);

    // masHelBrtDim
    masReadMeshDimensions(fileNames[2], "dim2", &masHelBrtDimMax[0]);
    masHelBrtDim = (float *)malloc(sizeof(float) * (int)(masHelBrtDimMax[0]));
    masReadMesh(fileNames[2], "dim2", &masHelBrtDim);

    // masHelBrrDim
    masReadMeshDimensions(fileNames[2], "dim1", &masHelBrrDimMax[0]);
    masHelBrrDim = (float *)malloc(sizeof(float) * (int)(masHelBrrDimMax[0]));
    masReadMesh(fileNames[2], "dim1", &masHelBrrDim);

    // grab the min and max for the r
    rTemp = masHelBrrDim[0];
    if (rTemp > rMin)
      rMin = rTemp;

    rTemp = masHelBrrDim[masHelBrrDimMax[0] - 1];
    if (rTemp < rMax)
      rMax = rTemp;

    // masHelVppDim
    masReadMeshDimensions(fileNames[3], "dim3", &masHelVppDimMax[0]);
    masHelVppDim = (float *)malloc(sizeof(float) * (int)(masHelVppDimMax[0]));
    masReadMesh(fileNames[3], "dim3", &masHelVppDim);

    // masHelVptDim
    masReadMeshDimensions(fileNames[3], "dim2", &masHelVptDimMax[0]);
    masHelVptDim = (float *)malloc(sizeof(float) * (int)(masHelVptDimMax[0]));
    masReadMesh(fileNames[3], "dim2", &masHelVptDim);

    // masHelVprDim
    masReadMeshDimensions(fileNames[3], "dim1", &masHelVprDimMax[0]);
    masHelVprDim = (float *)malloc(sizeof(float) * (int)(masHelVprDimMax[0]));
    masReadMesh(fileNames[3], "dim1", &masHelVprDim);

    // grab the min and max for the r
    rTemp = masHelVprDim[0];
    if (rTemp > rMin)
      rMin = rTemp;

    rTemp = masHelVprDim[masHelVprDimMax[0] - 1];
    if (rTemp < rMax)
      rMax = rTemp;

    // masHelVtpDim
    masReadMeshDimensions(fileNames[4], "dim3", &masHelVtpDimMax[0]);
    masHelVtpDim = (float *)malloc(sizeof(float) * (int)(masHelVtpDimMax[0]));
    masReadMesh(fileNames[4], "dim3", &masHelVtpDim);

    // masHelVttDim
    masReadMeshDimensions(fileNames[4], "dim2", &masHelVttDimMax[0]);
    masHelVttDim = (float *)malloc(sizeof(float) * (int)(masHelVttDimMax[0]));
    masReadMesh(fileNames[4], "dim2", &masHelVttDim);

    // masHelVtrDim
    masReadMeshDimensions(fileNames[4], "dim1", &masHelVtrDimMax[0]);
    masHelVtrDim = (float *)malloc(sizeof(float) * (int)(masHelVtrDimMax[0]));
    masReadMesh(fileNames[4], "dim1", &masHelVtrDim);

    // grab the min and max for the r
    rTemp = masHelVtrDim[0];
    if (rTemp > rMin)
      rMin = rTemp;

    rTemp = masHelVtrDim[masHelVtrDimMax[0] - 1];
    if (rTemp < rMax)
      rMax = rTemp;

    // masHelVrpDim
    masReadMeshDimensions(fileNames[5], "dim3", &masHelVrpDimMax[0]);
    masHelVrpDim = (float *)malloc(sizeof(float) * (int)(masHelVrpDimMax[0]));
    masReadMesh(fileNames[5], "dim3", &masHelVrpDim);

    // masHelVrtDim
    masReadMeshDimensions(fileNames[5], "dim2", &masHelVrtDimMax[0]);
    masHelVrtDim = (float *)malloc(sizeof(float) * (int)(masHelVrtDimMax[0]));
    masReadMesh(fileNames[5], "dim2", &masHelVrtDim);

    // masHelVrrDim
    masReadMeshDimensions(fileNames[5], "dim1", &masHelVrrDimMax[0]);
    masHelVrrDim = (float *)malloc(sizeof(float) * (int)(masHelVrrDimMax[0]));
    masReadMesh(fileNames[5], "dim1", &masHelVrrDim);

    // grab the min and max for the r
    rTemp = masHelVrrDim[0];
    if (rTemp > rMin)
      rMin = rTemp;

    rTemp = masHelVrrDim[masHelVrrDimMax[0] - 1];
    if (rTemp < rMax)
      rMax = rTemp;

    // masHelDpDim
    masReadMeshDimensions(fileNames[6], "dim3", &masHelDpDimMax[0]);
    masHelDpDim = (float *)malloc(sizeof(float) * (int)(masHelDpDimMax[0]));
    masReadMesh(fileNames[6], "dim3", &masHelDpDim);

    // masHelDtDim
    masReadMeshDimensions(fileNames[6], "dim2", &masHelDtDimMax[0]);
    masHelDtDim = (float *)malloc(sizeof(float) * (int)(masHelDtDimMax[0]));
    masReadMesh(fileNames[6], "dim2", &masHelDtDim);

    // masHelDrDim
    masReadMeshDimensions(fileNames[6], "dim1", &masHelDrDimMax[0]);
    masHelDrDim = (float *)malloc(sizeof(float) * (int)(masHelDrDimMax[0]));
    masReadMesh(fileNames[6], "dim1", &masHelDrDim);

    // grab the min and max for the r
    rTemp = masHelDrDim[0];
    if (rTemp > rMin)
      rMin = rTemp;

    rTemp = masHelDrDim[masHelDrDimMax[0] - 1];
    if (rTemp < rMax)
      rMax = rTemp;

  // set rScale, and masHelRadialMin/Max
  config.masHelRadialMin = rMin;
  config.masHelRadialMax = rMax;

  timer_mas_io = timer_mas_io + (MPI_Wtime() - timer_tmp);

}/*-------- END masHelReadFieldIndex()  ----------------------------*/
/*------------------------------------------------------------------*/

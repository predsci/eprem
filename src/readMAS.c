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
Index_t masResolutionChangeFlag;

int32_t masDimMin_0[1] = {0};
int32_t masDimMin_1[1] = {0};

int32_t masBppDimMax_0[1];
int32_t masBptDimMax_0[1];
int32_t masBprDimMax_0[1];
int32_t masHelBppDimMax_0[1];
int32_t masHelBptDimMax_0[1];
int32_t masHelBprDimMax_0[1];

int32_t masBppDimMax_1[1];
int32_t masBptDimMax_1[1];
int32_t masBprDimMax_1[1];
int32_t masHelBppDimMax_1[1];
int32_t masHelBptDimMax_1[1];
int32_t masHelBprDimMax_1[1];

int32_t masBtpDimMax_0[1];
int32_t masBttDimMax_0[1];
int32_t masBtrDimMax_0[1];
int32_t masHelBtpDimMax_0[1];
int32_t masHelBttDimMax_0[1];
int32_t masHelBtrDimMax_0[1];

int32_t masBtpDimMax_1[1];
int32_t masBttDimMax_1[1];
int32_t masBtrDimMax_1[1];
int32_t masHelBtpDimMax_1[1];
int32_t masHelBttDimMax_1[1];
int32_t masHelBtrDimMax_1[1];

int32_t masBrpDimMax_0[1];
int32_t masBrtDimMax_0[1];
int32_t masBrrDimMax_0[1];
int32_t masHelBrpDimMax_0[1];
int32_t masHelBrtDimMax_0[1];
int32_t masHelBrrDimMax_0[1];

int32_t masBrpDimMax_1[1];
int32_t masBrtDimMax_1[1];
int32_t masBrrDimMax_1[1];
int32_t masHelBrpDimMax_1[1];
int32_t masHelBrtDimMax_1[1];
int32_t masHelBrrDimMax_1[1];

int32_t masVppDimMax_0[1];
int32_t masVptDimMax_0[1];
int32_t masVprDimMax_0[1];
int32_t masHelVppDimMax_0[1];
int32_t masHelVptDimMax_0[1];
int32_t masHelVprDimMax_0[1];

int32_t masVppDimMax_1[1];
int32_t masVptDimMax_1[1];
int32_t masVprDimMax_1[1];
int32_t masHelVppDimMax_1[1];
int32_t masHelVptDimMax_1[1];
int32_t masHelVprDimMax_1[1];

int32_t masVtpDimMax_0[1];
int32_t masVttDimMax_0[1];
int32_t masVtrDimMax_0[1];
int32_t masHelVtpDimMax_0[1];
int32_t masHelVttDimMax_0[1];
int32_t masHelVtrDimMax_0[1];

int32_t masVtpDimMax_1[1];
int32_t masVttDimMax_1[1];
int32_t masVtrDimMax_1[1];
int32_t masHelVtpDimMax_1[1];
int32_t masHelVttDimMax_1[1];
int32_t masHelVtrDimMax_1[1];

int32_t masVrpDimMax_0[1];
int32_t masVrtDimMax_0[1];
int32_t masVrrDimMax_0[1];
int32_t masHelVrpDimMax_0[1];
int32_t masHelVrtDimMax_0[1];
int32_t masHelVrrDimMax_0[1];

int32_t masVrpDimMax_1[1];
int32_t masVrtDimMax_1[1];
int32_t masVrrDimMax_1[1];
int32_t masHelVrpDimMax_1[1];
int32_t masHelVrtDimMax_1[1];
int32_t masHelVrrDimMax_1[1];

int32_t masDpDimMax_0[1];
int32_t masDtDimMax_0[1];
int32_t masDrDimMax_0[1];
int32_t masHelDpDimMax_0[1];
int32_t masHelDtDimMax_0[1];
int32_t masHelDrDimMax_0[1];

int32_t masDpDimMax_1[1];
int32_t masDtDimMax_1[1];
int32_t masDrDimMax_1[1];
int32_t masHelDpDimMax_1[1];
int32_t masHelDtDimMax_1[1];
int32_t masHelDrDimMax_1[1];

// Mesh storage state0
float *masBppDim_0;
float *masBptDim_0;
float *masBprDim_0;
float *masHelBppDim_0;
float *masHelBptDim_0;
float *masHelBprDim_0;

float *masBtpDim_0;
float *masBttDim_0;
float *masBtrDim_0;
float *masHelBtpDim_0;
float *masHelBttDim_0;
float *masHelBtrDim_0;

float *masBrpDim_0;
float *masBrtDim_0;
float *masBrrDim_0;
float *masHelBrpDim_0;
float *masHelBrtDim_0;
float *masHelBrrDim_0;

float *masVppDim_0;
float *masVptDim_0;
float *masVprDim_0;
float *masHelVppDim_0;
float *masHelVptDim_0;
float *masHelVprDim_0;

float *masVtpDim_0;
float *masVttDim_0;
float *masVtrDim_0;
float *masHelVtpDim_0;
float *masHelVttDim_0;
float *masHelVtrDim_0;

float *masVrpDim_0;
float *masVrtDim_0;
float *masVrrDim_0;
float *masHelVrpDim_0;
float *masHelVrtDim_0;
float *masHelVrrDim_0;

float *masDpDim_0;
float *masDtDim_0;
float *masDrDim_0;
float *masHelDpDim_0;
float *masHelDtDim_0;
float *masHelDrDim_0;

// Mesh storage state1
float *masBppDim_1;
float *masBptDim_1;
float *masBprDim_1;
float *masHelBppDim_1;
float *masHelBptDim_1;
float *masHelBprDim_1;

float *masBtpDim_1;
float *masBttDim_1;
float *masBtrDim_1;
float *masHelBtpDim_1;
float *masHelBttDim_1;
float *masHelBtrDim_1;

float *masBrpDim_1;
float *masBrtDim_1;
float *masBrrDim_1;
float *masHelBrpDim_1;
float *masHelBrtDim_1;
float *masHelBrrDim_1;

float *masVppDim_1;
float *masVptDim_1;
float *masVprDim_1;
float *masHelVppDim_1;
float *masHelVptDim_1;
float *masHelVprDim_1;

float *masVtpDim_1;
float *masVttDim_1;
float *masVtrDim_1;
float *masHelVtpDim_1;
float *masHelVttDim_1;
float *masHelVtrDim_1;

float *masVrpDim_1;
float *masVrtDim_1;
float *masVrrDim_1;
float *masHelVrpDim_1;
float *masHelVrtDim_1;
float *masHelVrrDim_1;

float *masDpDim_1;
float *masDtDim_1;
float *masDrDim_1;
float *masHelDpDim_1;
float *masHelDtDim_1;
float *masHelDrDim_1;

// Data storage
float *masBp_0;
float *masBt_0;
float *masBr_0;
float *masVp_0;
float *masVt_0;
float *masVr_0;
float *masD_0;
float *masHelBp_0;
float *masHelBt_0;
float *masHelBr_0;
float *masHelVp_0;
float *masHelVt_0;
float *masHelVr_0;
float *masHelD_0;

float *masBp_1;
float *masBt_1;
float *masBr_1;
float *masVp_1;
float *masVt_1;
float *masVr_1;
float *masD_1;
float *masHelBp_1;
float *masHelBt_1;
float *masHelBr_1;
float *masHelVp_1;
float *masHelVt_1;
float *masHelVr_1;
float *masHelD_1;

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
    // masDatafile_type(); // This isnt working right - see other note.

    strncpy(file_extension,".h5",strlen(".h5")+1);

    masFileIndex0 = 0;
    masFileIndex1 = 0; // Check the initial state1 index
    // Initialize the dimensions and mesh for state0 and state1
    masReadFieldIndex("state0", masFileIndex0, masMallocFlag); //
    masReadFieldIndex("state1", masFileIndex1, masMallocFlag); //
    // heliospheric coupling
    if (config.masHelCouple > 0)
    {
      masHelReadFieldIndex("state0", masFileIndex0, masMallocFlag);
      masHelReadFieldIndex("state1", masFileIndex1, masMallocFlag); //
    }
    allocateMPIWindows("state0");
    allocateMPIWindows("state1");
    allocateHelMPIWindows("state0");
    allocateHelMPIWindows("state1");
    masMallocFlag = 1;

  }


}/*-------- END masFetchFileList()  --------------------------------*/
/*------------------------------------------------------------------*/

/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/*--*/     void                                                 /*--*/
/*--*/     cleanupMPIWindows(char *state)                       /*--*/
/*--                                                              --*/
/*--  This function cleans up MPI windows.                        --*/
/*------------------------------------------------------------------*/
{/*-----------------------------------------------------------------*/

   if (config.masCouple > 0) {
    if (strcmp(state, "state0") == 0)
    {
      MPI_Win_unlock_all(masBp_0_win);
      MPI_Win_unlock_all(masBt_0_win);
      MPI_Win_unlock_all(masBr_0_win);
      MPI_Win_unlock_all(masVp_0_win);
      MPI_Win_unlock_all(masVt_0_win);
      MPI_Win_unlock_all(masVr_0_win);
      MPI_Win_unlock_all(masD_0_win);
      MPI_Win_free(&masBp_0_win);
      MPI_Win_free(&masBt_0_win);
      MPI_Win_free(&masBr_0_win);
      MPI_Win_free(&masVp_0_win);
      MPI_Win_free(&masVt_0_win);
      MPI_Win_free(&masVr_0_win);
      MPI_Win_free(&masD_0_win);
    }
    else if (strcmp(state, "state1") == 0)
    {
      MPI_Win_unlock_all(masBp_1_win);
      MPI_Win_unlock_all(masBt_1_win);
      MPI_Win_unlock_all(masBr_1_win);
      MPI_Win_unlock_all(masVp_1_win);
      MPI_Win_unlock_all(masVt_1_win);
      MPI_Win_unlock_all(masVr_1_win);
      MPI_Win_unlock_all(masD_1_win);
      MPI_Win_free(&masBp_1_win);
      MPI_Win_free(&masBt_1_win);
      MPI_Win_free(&masBr_1_win);
      MPI_Win_free(&masVp_1_win);
      MPI_Win_free(&masVt_1_win);
      MPI_Win_free(&masVr_1_win);
      MPI_Win_free(&masD_1_win);
    }
  }

} /*-------- END cleanupMPIWindows()  --------------------------------*/
  /*------------------------------------------------------------------*/

/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/*--*/ void                       /*--*/
                                  /*--*/
cleanupHelMPIWindows(char *state) /*--*/
/*--                                                              --*/
/*--  This function cleans up Helio MPI windows.                  --*/
/*------------------------------------------------------------------*/
{ /*----------------------------------------------------------------*/
  // heliosphericoupling
  if (config.masHelCouple > 0)
  {
    if (strcmp(state, "state0") == 0)
    {
      MPI_Win_unlock_all(masHelBp_0_win);
      MPI_Win_unlock_all(masHelBt_0_win);
      MPI_Win_unlock_all(masHelBr_0_win);
      MPI_Win_unlock_all(masHelVp_0_win);
      MPI_Win_unlock_all(masHelVt_0_win);
      MPI_Win_unlock_all(masHelVr_0_win);
      MPI_Win_unlock_all(masHelD_0_win);

      MPI_Win_free(&masHelBp_0_win);
      MPI_Win_free(&masHelBt_0_win);
      MPI_Win_free(&masHelBr_0_win);
      MPI_Win_free(&masHelVp_0_win);
      MPI_Win_free(&masHelVt_0_win);
      MPI_Win_free(&masHelVr_0_win);
      MPI_Win_free(&masHelD_0_win);
    }
    else if (strcmp(state, "state1") == 0)
    {
      MPI_Win_unlock_all(masHelBp_1_win);
      MPI_Win_unlock_all(masHelBt_1_win);
      MPI_Win_unlock_all(masHelBr_1_win);
      MPI_Win_unlock_all(masHelVp_1_win);
      MPI_Win_unlock_all(masHelVt_1_win);
      MPI_Win_unlock_all(masHelVr_1_win);
      MPI_Win_unlock_all(masHelD_1_win);

      MPI_Win_free(&masHelBp_1_win);
      MPI_Win_free(&masHelBt_1_win);
      MPI_Win_free(&masHelBr_1_win);
      MPI_Win_free(&masHelVp_1_win);
      MPI_Win_free(&masHelVt_1_win);
      MPI_Win_free(&masHelVr_1_win);
      MPI_Win_free(&masHelD_1_win);
    }
  }

}/*-------- END cleanupHelMPIWindows()  ----------------------------*/
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
  int32_t masBppDimMax_test;

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
      if (masBppDimMax_0[0] == masBppDimMax_1[0]){ // state0 and previous state1 have the same dimension; just copy
        if (mpi_rank_shared==0) {
          N=(int)masBprDimMax_1[0]*(int)masBptDimMax_1[0]*(int)masBppDimMax_1[0];
          memcpy(&masBp_0[0],&masBp_1[0],N*sizeof(float));
          N=(int)masBtrDimMax_1[0]*(int)masBttDimMax_1[0]*(int)masBtpDimMax_1[0];
          memcpy(&masBt_0[0],&masBt_1[0],N*sizeof(float));
          N=(int)masBrrDimMax_1[0]*(int)masBrtDimMax_1[0]*(int)masBrpDimMax_1[0];
          memcpy(&masBr_0[0],&masBr_1[0],N*sizeof(float));
          N=(int)masVprDimMax_1[0]*(int)masVptDimMax_1[0]*(int)masVppDimMax_1[0];
          memcpy(&masVp_0[0],&masVp_1[0],N*sizeof(float));
          N=(int)masVtrDimMax_1[0]*(int)masVttDimMax_1[0]*(int)masVtpDimMax_1[0];
          memcpy(&masVt_0[0],&masVt_1[0],N*sizeof(float));
          N=(int)masVrrDimMax_1[0]*(int)masVrtDimMax_1[0]*(int)masVrpDimMax_1[0];
          memcpy(&masVr_0[0],&masVr_1[0],N*sizeof(float));
          N=(int)masDrDimMax_1[0] * (int)masDtDimMax_1[0]*(int)masDpDimMax_1[0];
          memcpy(&masD_0[0],&masD_1[0],N*sizeof(float));
        }
      }
      else { // Change in resolution detected; Reallocate state0 and copy from previous state1
        //printf("Resolution change detected in file id: %d; Reallocating state0 \n", masFileIndex0 + 1);
        // Deallocate and read in the dimensions and mesh for state0
        masReadFieldIndex("state0", masFileIndex0, masMallocFlag); //
        cleanupMPIWindows("state0");
        allocateMPIWindows("state0");
        // Now can copy from state1
        if (mpi_rank_shared == 0){
          N = (int)masBprDimMax_1[0] * (int)masBptDimMax_1[0] * (int)masBppDimMax_1[0];
          for (i = 0; i < N; i++){
            masBp_0[i] = masBp_1[i];
          }
          N = (int)masBtrDimMax_1[0] * (int)masBttDimMax_1[0] * (int)masBtpDimMax_1[0];
          for (i = 0; i < N; i++){
            masBt_0[i] = masBt_1[i];
          }
          N = (int)masBrrDimMax_1[0] * (int)masBrtDimMax_1[0] * (int)masBrpDimMax_1[0];
          for (i = 0; i < N; i++){
            masBr_0[i] = masBr_1[i];
          }
          N = (int)masVprDimMax_1[0] * (int)masVptDimMax_1[0] * (int)masVppDimMax_1[0];
          for (i = 0; i < N; i++){
            masVp_0[i] = masVp_1[i];
          }
          N = (int)masVtrDimMax_1[0] * (int)masVttDimMax_1[0] * (int)masVtpDimMax_1[0];
          for (i = 0; i < N; i++){
            masVt_0[i] = masVt_1[i];
          }
          N = (int)masVrrDimMax_1[0] * (int)masVrtDimMax_1[0] * (int)masVrpDimMax_1[0];
          for (i = 0; i < N; i++){
            masVr_0[i] = masVr_1[i];
          }
          N = (int)masDrDimMax_1[0] * (int)masDtDimMax_1[0] * (int)masDpDimMax_1[0];
          for (i = 0; i < N; i++){
            masD_0[i] = masD_1[i];
          }
        }
      }
    } else { // New data to be read into state0
      masBppDimMax_test = masReadBppDimTest(masFileIndex0); // Check for changes in grid resolution

      if (masBppDimMax_0[0] == masBppDimMax_test){
        if (mpi_rank_shared == 0){
          masReadData(masFileIndex0,
                  &masBp_0, &masBt_0, &masBr_0,
                  &masVp_0, &masVt_0, &masVr_0,
                  &masD_0);
        }
      }
      else { // Change in resolution detected
        printf("Resolution change detected in file id: %d; Reallocating state0 \n", masFileIndex0 + 1);
        // Deallocate and read in the dimensions and mesh for state0
        masReadFieldIndex("state0", masFileIndex0, masMallocFlag); //
        cleanupMPIWindows("state0");
        allocateMPIWindows("state0");
        if (mpi_rank_shared==0){
          masReadData(masFileIndex0,
                  &masBp_0, &masBt_0, &masBr_0,
                  &masVp_0, &masVt_0, &masVr_0,
                  &masD_0);
        }
      }
    }
    masFileIndex_loaded0 = masFileIndex0;
    need_sync_0 = 1;
  } else {
    need_sync_0 = 0;
  }

  if (masFileIndex1 != masFileIndex_loaded1){ // If state1 needs to be updated
    MPI_Barrier(comm_shared);
    masBppDimMax_test = masReadBppDimTest(masFileIndex1);  // Check if we need to read in the mesh again
    if (masBppDimMax_test == masBppDimMax_1[0])               // No change in grid resolution required
    {
      if (mpi_rank_shared == 0){
        masReadData(masFileIndex1,
                    &masBp_1, &masBt_1, &masBr_1,
                    &masVp_1, &masVt_1, &masVr_1,
                    &masD_1);
      }
    } else {
       // Need to re-read the state1 mesh
      printf("Resolution change detected in file id: %d; Reallocating state1 \n", masFileIndex1 + 1);
      // Deallocate and read in the dimensions and mesh for state0
      masReadFieldIndex("state1", masFileIndex1, masMallocFlag); //
      cleanupMPIWindows("state1");
      allocateMPIWindows("state1");
      if(mpi_rank_shared==0){
        masReadData(masFileIndex1,
                    &masBp_1, &masBt_1, &masBr_1,
                    &masVp_1, &masVt_1, &masVr_1,
                    &masD_1);
      }
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
          N=(int)masHelBprDimMax_0[0]*(int)masHelBptDimMax_0[0]*(int)masHelBppDimMax_0[0];
          memcpy(&masHelBp_0[0],&masHelBp_1[0],N*sizeof(float));
          N=(int)masHelBtrDimMax_0[0]*(int)masHelBttDimMax_0[0]*(int)masHelBtpDimMax_0[0];
          memcpy(&masHelBt_0[0],&masHelBt_1[0],N*sizeof(float));
          N=(int)masHelBrrDimMax_0[0]*(int)masHelBrtDimMax_0[0]*(int)masHelBrpDimMax_0[0];
          memcpy(&masHelBr_0[0],&masHelBr_1[0],N*sizeof(float));
          N=(int)masHelVprDimMax_0[0]*(int)masHelVptDimMax_0[0]*(int)masHelVppDimMax_0[0];
          memcpy(&masHelVp_0[0],&masHelVp_1[0],N*sizeof(float));
          N=(int)masHelVtrDimMax_0[0]*(int)masHelVttDimMax_0[0]*(int)masHelVtpDimMax_0[0];
          memcpy(&masHelVt_0[0],&masHelVt_1[0],N*sizeof(float));
          N=(int)masHelVrrDimMax_0[0]*(int)masHelVrtDimMax_0[0]*(int)masHelVrpDimMax_0[0];
          memcpy(&masHelVr_0[0],&masHelVr_1[0],N*sizeof(float));
          N=(int)masHelDrDimMax_0[0]*(int)masHelDtDimMax_0[0]*(int)masHelDpDimMax_0[0];
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
/*--*/     masReadFieldIndex(char *state, Index_t fileIndex, Index_t MallocFlag) /*--*/
/*--                                                              --*/
/*------------------------------------------------------------------*/
{/*-----------------------------------------------------------------*/

  char fileNames[7][MAX_STRING_SIZE];

  Scalar_t rMin, rMax, rTemp;

  double timer_tmp=0;

  timer_tmp = MPI_Wtime();

  // Deallocate if needed
  if (MallocFlag == 1)
  {
    if (strcmp(state, "state0") == 0)
    {
      free(masBppDim_0);
      free(masBptDim_0);
      free(masBprDim_0);
      free(masBtpDim_0);
      free(masBttDim_0);
      free(masBtrDim_0);
      free(masBrpDim_0);
      free(masBrtDim_0);
      free(masBrrDim_0);

      free(masVppDim_0);
      free(masVptDim_0);
      free(masVprDim_0);
      free(masVtpDim_0);
      free(masVttDim_0);
      free(masVtrDim_0);
      free(masVrpDim_0);
      free(masVrtDim_0);
      free(masVrrDim_0);

      free(masDpDim_0);
      free(masDtDim_0);
      free(masDrDim_0);
    }
    else if (strcmp(state, "state1") == 0)
    {
      free(masBppDim_1);
      free(masBptDim_1);
      free(masBprDim_1);
      free(masBtpDim_1);
      free(masBttDim_1);
      free(masBtrDim_1);
      free(masBrpDim_1);
      free(masBrtDim_1);
      free(masBrrDim_1);

      free(masVppDim_1);
      free(masVptDim_1);
      free(masVprDim_1);
      free(masVtpDim_1);
      free(masVttDim_1);
      free(masVtrDim_1);
      free(masVrpDim_1);
      free(masVrtDim_1);
      free(masVrrDim_1);

      free(masDpDim_1);
      free(masDtDim_1);
      free(masDrDim_1);
    }
    else
    {
      printf("Unknown state \n");
      exit(0);
    }
  }

  // Read in dimensions --> Allocate the mesh storage --> Read the mesh

  if (config.masDigits == 3)
  {
    sprintf(fileNames[0], "%sbp%03d%s", config.masDirectory, fileIndex + 1, file_extension);
    sprintf(fileNames[1], "%sbt%03d%s", config.masDirectory, fileIndex + 1, file_extension);
    sprintf(fileNames[2], "%sbr%03d%s", config.masDirectory, fileIndex + 1, file_extension);
    sprintf(fileNames[3], "%svp%03d%s", config.masDirectory, fileIndex + 1, file_extension);
    sprintf(fileNames[4], "%svt%03d%s", config.masDirectory, fileIndex + 1, file_extension);
    sprintf(fileNames[5], "%svr%03d%s", config.masDirectory, fileIndex + 1, file_extension);
    sprintf(fileNames[6], "%srho%03d%s", config.masDirectory, fileIndex + 1, file_extension);
  }
  else
  {
    sprintf(fileNames[0], "%sbp%06d%s", config.masDirectory, fileIndex + 1, file_extension);
    sprintf(fileNames[1], "%sbt%06d%s", config.masDirectory, fileIndex + 1, file_extension);
    sprintf(fileNames[2], "%sbr%06d%s", config.masDirectory, fileIndex + 1, file_extension);
    sprintf(fileNames[3], "%svp%06d%s", config.masDirectory, fileIndex + 1, file_extension);
    sprintf(fileNames[4], "%svt%06d%s", config.masDirectory, fileIndex + 1, file_extension);
    sprintf(fileNames[5], "%svr%06d%s", config.masDirectory, fileIndex + 1, file_extension);
    sprintf(fileNames[6], "%srho%06d%s", config.masDirectory, fileIndex + 1, file_extension);
  }

  //

  if (strcmp(state, "state0") == 0)
  {
    // masBpp
    masReadMeshDimensions(fileNames[0], "dim3", &masBppDimMax_0[0]);
    masBppDim_0 = (float *)malloc(sizeof(float) * (int)(masBppDimMax_0[0]));
    masReadMesh(fileNames[0], "dim3", &masBppDim_0);

    // masBptDim
    masReadMeshDimensions(fileNames[0], "dim2", &masBptDimMax_0[0]);
    masBptDim_0 = (float *)malloc(sizeof(float) * (int)(masBptDimMax_0[0]));
    masReadMesh(fileNames[0], "dim2", &masBptDim_0);

    // masBprDim
    masReadMeshDimensions(fileNames[0], "dim1", &masBprDimMax_0[0]);
    masBprDim_0 = (float *)malloc(sizeof(float) * (int)(masBprDimMax_0[0]));
    masReadMesh(fileNames[0], "dim1", &masBprDim_0);

    // grab the min and max for the r
    rMin = masBprDim_0[0];
    rMax = masBprDim_0[masBprDimMax_0[0] - 1];

    // masBtpDim
    masReadMeshDimensions(fileNames[1], "dim3", &masBtpDimMax_0[0]);
    masBtpDim_0 = (float *)malloc(sizeof(float) * (int)(masBtpDimMax_0[0]));
    masReadMesh(fileNames[1], "dim3", &masBtpDim_0);

    // masBttDim
    masReadMeshDimensions(fileNames[1], "dim2", &masBttDimMax_0[0]);
    masBttDim_0 = (float *)malloc(sizeof(float) * (int)(masBttDimMax_0[0]));
    masReadMesh(fileNames[1], "dim2", &masBttDim_0);

    // masBtrDim
    masReadMeshDimensions(fileNames[1], "dim1", &masBtrDimMax_0[0]);
    masBtrDim_0 = (float *)malloc(sizeof(float) * (int)(masBtrDimMax_0[0]));
    masReadMesh(fileNames[1], "dim1", &masBtrDim_0);

    // grab the min and max for the r
    rTemp = masBtrDim_0[0];
    if (rTemp > rMin)
      rMin = rTemp;

    rTemp = masBtrDim_0[masBtrDimMax_0[0] - 1];
    if (rTemp < rMax)
      rMax = rTemp;

    // masBrpDim
    masReadMeshDimensions(fileNames[2], "dim3", &masBrpDimMax_0[0]);
    masBrpDim_0 = (float *)malloc(sizeof(float) * (int)(masBrpDimMax_0[0]));
    masReadMesh(fileNames[2], "dim3", &masBrpDim_0);

    // masBrtDim
    masReadMeshDimensions(fileNames[2], "dim2", &masBrtDimMax_0[0]);
    masBrtDim_0 = (float *)malloc(sizeof(float) * (int)(masBrtDimMax_0[0]));
    masReadMesh(fileNames[2], "dim2", &masBrtDim_0);

    // masBrrDim
    masReadMeshDimensions(fileNames[2], "dim1", &masBrrDimMax_0[0]);
    masBrrDim_0 = (float *)malloc(sizeof(float) * (int)(masBrrDimMax_0[0]));
    masReadMesh(fileNames[2], "dim1", &masBrrDim_0);

    // grab the min and max for the r
    rTemp = masBrrDim_0[0];
    if (rTemp > rMin)
      rMin = rTemp;

    rTemp = masBrrDim_0[masBrrDimMax_0[0] - 1];
    if (rTemp < rMax)
      rMax = rTemp;

    // masVppDim
    masReadMeshDimensions(fileNames[3], "dim3", &masVppDimMax_0[0]);
    masVppDim_0 = (float *)malloc(sizeof(float) * (int)(masVppDimMax_0[0]));
    masReadMesh(fileNames[3], "dim3", &masVppDim_0);

    // masVptDim
    masReadMeshDimensions(fileNames[3], "dim2", &masVptDimMax_0[0]);
    masVptDim_0 = (float *)malloc(sizeof(float) * (int)(masVptDimMax_0[0]));
    masReadMesh(fileNames[3], "dim2", &masVptDim_0);

    // masVprDim
    masReadMeshDimensions(fileNames[3], "dim1", &masVprDimMax_0[0]);
    masVprDim_0 = (float *)malloc(sizeof(float) * (int)(masVprDimMax_0[0]));
    masReadMesh(fileNames[3], "dim1", &masVprDim_0);

    // grab the min and max for the r
    rTemp = masVprDim_0[0];
    if (rTemp > rMin)
      rMin = rTemp;

    rTemp = masVprDim_0[masVprDimMax_0[0] - 1];
    if (rTemp < rMax)
      rMax = rTemp;

    // masVtpDim
    masReadMeshDimensions(fileNames[4], "dim3", &masVtpDimMax_0[0]);
    masVtpDim_0 = (float *)malloc(sizeof(float) * (int)(masVtpDimMax_0[0]));
    masReadMesh(fileNames[4], "dim3", &masVtpDim_0);

    // masVttDim
    masReadMeshDimensions(fileNames[4], "dim2", &masVttDimMax_0[0]);
    masVttDim_0 = (float *)malloc(sizeof(float) * (int)(masVttDimMax_0[0]));
    masReadMesh(fileNames[4], "dim2", &masVttDim_0);

    // masVtrDim
    masReadMeshDimensions(fileNames[4], "dim1", &masVtrDimMax_0[0]);
    masVtrDim_0 = (float *)malloc(sizeof(float) * (int)(masVtrDimMax_0[0]));
    masReadMesh(fileNames[4], "dim1", &masVtrDim_0);

    // grab the min and max for the r
    rTemp = masVtrDim_0[0];
    if (rTemp > rMin)
      rMin = rTemp;

    rTemp = masVtrDim_0[masVtrDimMax_0[0] - 1];
    if (rTemp < rMax)
      rMax = rTemp;

    // masVrpDim
    masReadMeshDimensions(fileNames[5], "dim3", &masVrpDimMax_0[0]);
    masVrpDim_0 = (float *)malloc(sizeof(float) * (int)(masVrpDimMax_0[0]));
    masReadMesh(fileNames[5], "dim3", &masVrpDim_0);

    // masVrtDim
    masReadMeshDimensions(fileNames[5], "dim2", &masVrtDimMax_0[0]);
    masVrtDim_0 = (float *)malloc(sizeof(float) * (int)(masVrtDimMax_0[0]));
    masReadMesh(fileNames[5], "dim2", &masVrtDim_0);

    // masVrrDim
    masReadMeshDimensions(fileNames[5], "dim1", &masVrrDimMax_0[0]);
    masVrrDim_0 = (float *)malloc(sizeof(float) * (int)(masVrrDimMax_0[0]));
    masReadMesh(fileNames[5], "dim1", &masVrrDim_0);

    // grab the min and max for the r
    rTemp = masVrrDim_0[0];
    if (rTemp > rMin)
      rMin = rTemp;

    rTemp = masVrrDim_0[masVrrDimMax_0[0] - 1];
    if (rTemp < rMax)
      rMax = rTemp;

    // masDpDim
    masReadMeshDimensions(fileNames[6], "dim3", &masDpDimMax_0[0]);
    masDpDim_0 = (float *)malloc(sizeof(float) * (int)(masDpDimMax_0[0]));
    masReadMesh(fileNames[6], "dim3", &masDpDim_0);

    // masDtDim
    masReadMeshDimensions(fileNames[6], "dim2", &masDtDimMax_0[0]);
    masDtDim_0 = (float *)malloc(sizeof(float) * (int)(masDtDimMax_0[0]));
    masReadMesh(fileNames[6], "dim2", &masDtDim_0);

    // masDrDim
    masReadMeshDimensions(fileNames[6], "dim1", &masDrDimMax_0[0]);
    masDrDim_0 = (float *)malloc(sizeof(float) * (int)(masDrDimMax_0[0]));
    masReadMesh(fileNames[6], "dim1", &masDrDim_0);

    // grab the min and max for the r
    rTemp = masDrDim_0[0];
    if (rTemp > rMin)
      rMin = rTemp;

    rTemp = masDrDim_0[masDrDimMax_0[0] - 1];
    if (rTemp < rMax)
      rMax = rTemp;
  }
  else if (strcmp(state, "state1") == 0)
  {
    // masBpp
    masReadMeshDimensions(fileNames[0], "dim3", &masBppDimMax_1[0]);
    masBppDim_1 = (float *)malloc(sizeof(float) * (int)(masBppDimMax_1[0]));
    masReadMesh(fileNames[0], "dim3", &masBppDim_1);

    // masBptDim
    masReadMeshDimensions(fileNames[0], "dim2", &masBptDimMax_1[0]);
    masBptDim_1 = (float *)malloc(sizeof(float) * (int)(masBptDimMax_1[0]));
    masReadMesh(fileNames[0], "dim2", &masBptDim_1);

    // masBprDim
    masReadMeshDimensions(fileNames[0], "dim1", &masBprDimMax_1[0]);
    masBprDim_1 = (float *)malloc(sizeof(float) * (int)(masBprDimMax_1[0]));
    masReadMesh(fileNames[0], "dim1", &masBprDim_1);

    // grab the min and max for the r
    rMin = masBprDim_1[0];
    rMax = masBprDim_1[masBprDimMax_1[0] - 1];

    // masBtpDim
    masReadMeshDimensions(fileNames[1], "dim3", &masBtpDimMax_1[0]);
    masBtpDim_1 = (float *)malloc(sizeof(float) * (int)(masBtpDimMax_1[0]));
    masReadMesh(fileNames[1], "dim3", &masBtpDim_1);

    // masBttDim
    masReadMeshDimensions(fileNames[1], "dim2", &masBttDimMax_1[0]);
    masBttDim_1 = (float *)malloc(sizeof(float) * (int)(masBttDimMax_1[0]));
    masReadMesh(fileNames[1], "dim2", &masBttDim_1);

    // masBtrDim
    masReadMeshDimensions(fileNames[1], "dim1", &masBtrDimMax_1[0]);
    masBtrDim_1 = (float *)malloc(sizeof(float) * (int)(masBtrDimMax_1[0]));
    masReadMesh(fileNames[1], "dim1", &masBtrDim_1);

    // grab the min and max for the r
    rTemp = masBtrDim_1[0];
    if (rTemp > rMin)
      rMin = rTemp;

    rTemp = masBtrDim_1[masBtrDimMax_1[0] - 1];
    if (rTemp < rMax)
      rMax = rTemp;

    // masBrpDim
    masReadMeshDimensions(fileNames[2], "dim3", &masBrpDimMax_1[0]);
    masBrpDim_1 = (float *)malloc(sizeof(float) * (int)(masBrpDimMax_1[0]));
    masReadMesh(fileNames[2], "dim3", &masBrpDim_1);

    // masBrtDim
    masReadMeshDimensions(fileNames[2], "dim2", &masBrtDimMax_1[0]);
    masBrtDim_1 = (float *)malloc(sizeof(float) * (int)(masBrtDimMax_1[0]));
    masReadMesh(fileNames[2], "dim2", &masBrtDim_1);

    // masBrrDim
    masReadMeshDimensions(fileNames[2], "dim1", &masBrrDimMax_1[0]);
    masBrrDim_1 = (float *)malloc(sizeof(float) * (int)(masBrrDimMax_1[0]));
    masReadMesh(fileNames[2], "dim1", &masBrrDim_1);

    // grab the min and max for the r
    rTemp = masBrrDim_1[0];
    if (rTemp > rMin)
      rMin = rTemp;

    rTemp = masBrrDim_1[masBrrDimMax_1[0] - 1];
    if (rTemp < rMax)
      rMax = rTemp;

    // masVppDim
    masReadMeshDimensions(fileNames[3], "dim3", &masVppDimMax_1[0]);
    masVppDim_1 = (float *)malloc(sizeof(float) * (int)(masVppDimMax_1[0]));
    masReadMesh(fileNames[3], "dim3", &masVppDim_1);

    // masVptDim
    masReadMeshDimensions(fileNames[3], "dim2", &masVptDimMax_1[0]);
    masVptDim_1 = (float *)malloc(sizeof(float) * (int)(masVptDimMax_1[0]));
    masReadMesh(fileNames[3], "dim2", &masVptDim_1);

    // masVprDim
    masReadMeshDimensions(fileNames[3], "dim1", &masVprDimMax_1[0]);
    masVprDim_1 = (float *)malloc(sizeof(float) * (int)(masVprDimMax_1[0]));
    masReadMesh(fileNames[3], "dim1", &masVprDim_1);

    // grab the min and max for the r
    rTemp = masVprDim_1[0];
    if (rTemp > rMin)
      rMin = rTemp;

    rTemp = masVprDim_1[masVprDimMax_1[0] - 1];
    if (rTemp < rMax)
      rMax = rTemp;

    // masVtpDim
    masReadMeshDimensions(fileNames[4], "dim3", &masVtpDimMax_1[0]);
    masVtpDim_1 = (float *)malloc(sizeof(float) * (int)(masVtpDimMax_1[0]));
    masReadMesh(fileNames[4], "dim3", &masVtpDim_1);

    // masVttDim
    masReadMeshDimensions(fileNames[4], "dim2", &masVttDimMax_1[0]);
    masVttDim_1 = (float *)malloc(sizeof(float) * (int)(masVttDimMax_1[0]));
    masReadMesh(fileNames[4], "dim2", &masVttDim_1);

    // masVtrDim
    masReadMeshDimensions(fileNames[4], "dim1", &masVtrDimMax_1[0]);
    masVtrDim_1 = (float *)malloc(sizeof(float) * (int)(masVtrDimMax_1[0]));
    masReadMesh(fileNames[4], "dim1", &masVtrDim_1);

    // grab the min and max for the r
    rTemp = masVtrDim_1[0];
    if (rTemp > rMin)
      rMin = rTemp;

    rTemp = masVtrDim_1[masVtrDimMax_1[0] - 1];
    if (rTemp < rMax)
      rMax = rTemp;

    // masVrpDim
    masReadMeshDimensions(fileNames[5], "dim3", &masVrpDimMax_1[0]);
    masVrpDim_1 = (float *)malloc(sizeof(float) * (int)(masVrpDimMax_1[0]));
    masReadMesh(fileNames[5], "dim3", &masVrpDim_1);

    // masVrtDim
    masReadMeshDimensions(fileNames[5], "dim2", &masVrtDimMax_1[0]);
    masVrtDim_1 = (float *)malloc(sizeof(float) * (int)(masVrtDimMax_1[0]));
    masReadMesh(fileNames[5], "dim2", &masVrtDim_1);

    // masVrrDim
    masReadMeshDimensions(fileNames[5], "dim1", &masVrrDimMax_1[0]);
    masVrrDim_1 = (float *)malloc(sizeof(float) * (int)(masVrrDimMax_1[0]));
    masReadMesh(fileNames[5], "dim1", &masVrrDim_1);

    // grab the min and max for the r
    rTemp = masVrrDim_1[0];
    if (rTemp > rMin)
      rMin = rTemp;

    rTemp = masVrrDim_1[masVrrDimMax_1[0] - 1];
    if (rTemp < rMax)
      rMax = rTemp;

    // masDpDim
    masReadMeshDimensions(fileNames[6], "dim3", &masDpDimMax_1[0]);
    masDpDim_1 = (float *)malloc(sizeof(float) * (int)(masDpDimMax_1[0]));
    masReadMesh(fileNames[6], "dim3", &masDpDim_1);

    // masDtDim
    masReadMeshDimensions(fileNames[6], "dim2", &masDtDimMax_1[0]);
    masDtDim_1 = (float *)malloc(sizeof(float) * (int)(masDtDimMax_1[0]));
    masReadMesh(fileNames[6], "dim2", &masDtDim_1);

    // masDrDim
    masReadMeshDimensions(fileNames[6], "dim1", &masDrDimMax_1[0]);
    masDrDim_1 = (float *)malloc(sizeof(float) * (int)(masDrDimMax_1[0]));
    masReadMesh(fileNames[6], "dim1", &masDrDim_1);

    // grab the min and max for the r
    rTemp = masDrDim_1[0];
    if (rTemp > rMin)
      rMin = rTemp;

    rTemp = masDrDim_1[masDrDimMax_1[0] - 1];
    if (rTemp < rMax)
      rMax = rTemp;
  }
  else
  {
    printf("Unknown state \n");
    exit(0);
  }

  // set rScale, and masRadialMin/Max (?) Cooper= config.rScale= RSAU
  config.rScale = rMin * RSAU;
  config.masRadialMin = rMin;
  config.masRadialMax = rMax;

  timer_mas_io = timer_mas_io + (MPI_Wtime() - timer_tmp);

} /*-------- END masReadFieldIndex()  -----------------------*/
  /*------------------------------------------------------------------*/

/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/*--*/ void                                                     /*--*/
                                                                /*--*/
masHelReadFieldIndex(char *state, Index_t fileIndex, Index_t MallocFlag) /*--*/
/*--                                                              --*/
/*------------------------------------------------------------------*/
{/*-----------------------------------------------------------------*/

  char fileNames[7][MAX_STRING_SIZE];
  Scalar_t rMin, rMax, rTemp;

  double timer_tmp = 0;

  timer_tmp = MPI_Wtime();

  // Deallocate if needed
  if (MallocFlag == 1)
  {
    if (strcmp(state,"state0")==0)
    {
      free(masHelBppDim_0);
      free(masHelBptDim_0);
      free(masHelBprDim_0);
      free(masHelBtpDim_0);
      free(masHelBttDim_0);
      free(masHelBtrDim_0);
      free(masHelBrpDim_0);
      free(masHelBrtDim_0);
      free(masHelBrrDim_0);

      free(masHelVppDim_0);
      free(masHelVptDim_0);
      free(masHelVprDim_0);
      free(masHelVtpDim_0);
      free(masHelVttDim_0);
      free(masHelVtrDim_0);
      free(masHelVrpDim_0);
      free(masHelVrtDim_0);
      free(masHelVrrDim_0);

      free(masHelDpDim_0);
      free(masHelDtDim_0);
      free(masHelDrDim_0);
    }
    else if (strcmp(state, "state1")==0)
    {
      free(masHelBppDim_1);
      free(masHelBptDim_1);
      free(masHelBprDim_1);
      free(masHelBtpDim_1);
      free(masHelBttDim_1);
      free(masHelBtrDim_1);
      free(masHelBrpDim_1);
      free(masHelBrtDim_1);
      free(masHelBrrDim_1);

      free(masHelVppDim_1);
      free(masHelVptDim_1);
      free(masHelVprDim_1);
      free(masHelVtpDim_1);
      free(masHelVttDim_1);
      free(masHelVtrDim_1);
      free(masHelVrpDim_1);
      free(masHelVrtDim_1);
      free(masHelVrrDim_1);

      free(masHelDpDim_1);
      free(masHelDtDim_1);
      free(masHelDrDim_1);
    }
    else
    {
      printf("Unknown state \n");
      exit(0);
    }
  }

  // Read in dimensions --> Allocate the mesh storage --> Read the mesh


  if (config.masHelDigits == 3)
  {
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

  if (strcmp(state, "state0")==0)
  {
    // masHelBpp
    masReadMeshDimensions(fileNames[0], "dim3", &masHelBppDimMax_0[0]);
    masHelBppDim_0 = (float *)malloc(sizeof(float) * (int)(masHelBppDimMax_0[0]));
    masReadMesh(fileNames[0], "dim3", &masHelBppDim_0);

    // masHelBptDim
    masReadMeshDimensions(fileNames[0], "dim2", &masHelBptDimMax_0[0]);
    masHelBptDim_0 = (float *)malloc(sizeof(float) * (int)(masHelBptDimMax_0[0]));
    masReadMesh(fileNames[0], "dim2", &masHelBptDim_0);

    // masHelBprDim
    masReadMeshDimensions(fileNames[0], "dim1", &masHelBprDimMax_0[0]);
    masHelBprDim_0 = (float *)malloc(sizeof(float) * (int)(masHelBprDimMax_0[0]));
    masReadMesh(fileNames[0], "dim1", &masHelBprDim_0);
    // grab the min and max for the r
    rMin = masHelBprDim_0[0];
    rMax = masHelBprDim_0[masHelBprDimMax_0[0] - 1];

    // masHelBtpDim
    masReadMeshDimensions(fileNames[1], "dim3", &masHelBtpDimMax_0[0]);
    masHelBtpDim_0 = (float *)malloc(sizeof(float) * (int)(masHelBtpDimMax_0[0]));
    masReadMesh(fileNames[1], "dim3", &masHelBtpDim_0);

    // masHelBttDim
    masReadMeshDimensions(fileNames[1], "dim2", &masHelBttDimMax_0[0]);
    masHelBttDim_0 = (float *)malloc(sizeof(float) * (int)(masHelBttDimMax_0[0]));
    masReadMesh(fileNames[1], "dim2", &masHelBttDim_0);

    // masHelBtrDim
    masReadMeshDimensions(fileNames[1], "dim1", &masHelBtrDimMax_0[0]);
    masHelBtrDim_0 = (float *)malloc(sizeof(float) * (int)(masHelBtrDimMax_0[0]));
    masReadMesh(fileNames[1], "dim1", &masHelBtrDim_0);

    // grab the min and max for the r
    rTemp = masHelBtrDim_0[0];
    if (rTemp > rMin)
      rMin = rTemp;

    rTemp = masHelBtrDim_0[masHelBtrDimMax_0[0] - 1];
    if (rTemp < rMax)
      rMax = rTemp;

    // masHelBrpDim
    masReadMeshDimensions(fileNames[2], "dim3", &masHelBrpDimMax_0[0]);
    masHelBrpDim_0 = (float *)malloc(sizeof(float) * (int)(masHelBrpDimMax_0[0]));
    masReadMesh(fileNames[2], "dim3", &masHelBrpDim_0);

    // masHelBrtDim
    masReadMeshDimensions(fileNames[2], "dim2", &masHelBrtDimMax_0[0]);
    masHelBrtDim_0 = (float *)malloc(sizeof(float) * (int)(masHelBrtDimMax_0[0]));
    masReadMesh(fileNames[2], "dim2", &masHelBrtDim_0);

    // masHelBrrDim
    masReadMeshDimensions(fileNames[2], "dim1", &masHelBrrDimMax_0[0]);
    masHelBrrDim_0 = (float *)malloc(sizeof(float) * (int)(masHelBrrDimMax_0[0]));
    masReadMesh(fileNames[2], "dim1", &masHelBrrDim_0);

    // grab the min and max for the r
    rTemp = masHelBrrDim_0[0];
    if (rTemp > rMin)
      rMin = rTemp;

    rTemp = masHelBrrDim_0[masHelBrrDimMax_0[0] - 1];
    if (rTemp < rMax)
      rMax = rTemp;

    // masHelVppDim
    masReadMeshDimensions(fileNames[3], "dim3", &masHelVppDimMax_0[0]);
    masHelVppDim_0 = (float *)malloc(sizeof(float) * (int)(masHelVppDimMax_0[0]));
    masReadMesh(fileNames[3], "dim3", &masHelVppDim_0);

    // masHelVptDim
    masReadMeshDimensions(fileNames[3], "dim2", &masHelVptDimMax_0[0]);
    masHelVptDim_0 = (float *)malloc(sizeof(float) * (int)(masHelVptDimMax_0[0]));
    masReadMesh(fileNames[3], "dim2", &masHelVptDim_0);

    // masHelVprDim
    masReadMeshDimensions(fileNames[3], "dim1", &masHelVprDimMax_0[0]);
    masHelVprDim_0 = (float *)malloc(sizeof(float) * (int)(masHelVprDimMax_0[0]));
    masReadMesh(fileNames[3], "dim1", &masHelVprDim_0);

    // grab the min and max for the r
    rTemp = masHelVprDim_0[0];
    if (rTemp > rMin)
      rMin = rTemp;

    rTemp = masHelVprDim_0[masHelVprDimMax_0[0] - 1];
    if (rTemp < rMax)
      rMax = rTemp;

    // masHelVtpDim
    masReadMeshDimensions(fileNames[4], "dim3", &masHelVtpDimMax_0[0]);
    masHelVtpDim_0 = (float *)malloc(sizeof(float) * (int)(masHelVtpDimMax_0[0]));
    masReadMesh(fileNames[4], "dim3", &masHelVtpDim_0);

    // masHelVttDim
    masReadMeshDimensions(fileNames[4], "dim2", &masHelVttDimMax_0[0]);
    masHelVttDim_0 = (float *)malloc(sizeof(float) * (int)(masHelVttDimMax_0[0]));
    masReadMesh(fileNames[4], "dim2", &masHelVttDim_0);

    // masHelVtrDim
    masReadMeshDimensions(fileNames[4], "dim1", &masHelVtrDimMax_0[0]);
    masHelVtrDim_0 = (float *)malloc(sizeof(float) * (int)(masHelVtrDimMax_0[0]));
    masReadMesh(fileNames[4], "dim1", &masHelVtrDim_0);

    // grab the min and max for the r
    rTemp = masHelVtrDim_0[0];
    if (rTemp > rMin)
      rMin = rTemp;

    rTemp = masHelVtrDim_0[masHelVtrDimMax_0[0] - 1];
    if (rTemp < rMax)
      rMax = rTemp;

    // masHelVrpDim
    masReadMeshDimensions(fileNames[5], "dim3", &masHelVrpDimMax_0[0]);
    masHelVrpDim_0 = (float *)malloc(sizeof(float) * (int)(masHelVrpDimMax_0[0]));
    masReadMesh(fileNames[5], "dim3", &masHelVrpDim_0);

    // masHelVrtDim
    masReadMeshDimensions(fileNames[5], "dim2", &masHelVrtDimMax_0[0]);
    masHelVrtDim_0 = (float *)malloc(sizeof(float) * (int)(masHelVrtDimMax_0[0]));
    masReadMesh(fileNames[5], "dim2", &masHelVrtDim_0);

    // masHelVrrDim
    masReadMeshDimensions(fileNames[5], "dim1", &masHelVrrDimMax_0[0]);
    masHelVrrDim_0 = (float *)malloc(sizeof(float) * (int)(masHelVrrDimMax_0[0]));
    masReadMesh(fileNames[5], "dim1", &masHelVrrDim_0);

    // grab the min and max for the r
    rTemp = masHelVrrDim_0[0];
    if (rTemp > rMin)
      rMin = rTemp;

    rTemp = masHelVrrDim_0[masHelVrrDimMax_0[0] - 1];
    if (rTemp < rMax)
      rMax = rTemp;

    // masHelDpDim
    masReadMeshDimensions(fileNames[6], "dim3", &masHelDpDimMax_0[0]);
    masHelDpDim_0 = (float *)malloc(sizeof(float) * (int)(masHelDpDimMax_0[0]));
    masReadMesh(fileNames[6], "dim3", &masHelDpDim_0);

    // masHelDtDim
    masReadMeshDimensions(fileNames[6], "dim2", &masHelDtDimMax_0[0]);
    masHelDtDim_0 = (float *)malloc(sizeof(float) * (int)(masHelDtDimMax_0[0]));
    masReadMesh(fileNames[6], "dim2", &masHelDtDim_0);

    // masHelDrDim
    masReadMeshDimensions(fileNames[6], "dim1", &masHelDrDimMax_0[0]);
    masHelDrDim_0 = (float *)malloc(sizeof(float) * (int)(masHelDrDimMax_0[0]));
    masReadMesh(fileNames[6], "dim1", &masHelDrDim_0);

    // grab the min and max for the r
    rTemp = masHelDrDim_0[0];
    if (rTemp > rMin)
      rMin = rTemp;

    rTemp = masHelDrDim_0[masHelDrDimMax_0[0] - 1];
    if (rTemp < rMax)
      rMax = rTemp;
  }
  else if (strcmp(state, "state1")==0)
  {
    // masHelBpp
    masReadMeshDimensions(fileNames[0], "dim3", &masHelBppDimMax_1[0]);
    masHelBppDim_1 = (float *)malloc(sizeof(float) * (int)(masHelBppDimMax_1[0]));
    masReadMesh(fileNames[0], "dim3", &masHelBppDim_1);

    // masHelBptDim
    masReadMeshDimensions(fileNames[0], "dim2", &masHelBptDimMax_1[0]);
    masHelBptDim_1 = (float *)malloc(sizeof(float) * (int)(masHelBptDimMax_1[0]));
    masReadMesh(fileNames[0], "dim2", &masHelBptDim_1);

    // masHelBprDim
    masReadMeshDimensions(fileNames[0], "dim1", &masHelBprDimMax_1[0]);
    masHelBprDim_1 = (float *)malloc(sizeof(float) * (int)(masHelBprDimMax_1[0]));
    masReadMesh(fileNames[0], "dim1", &masHelBprDim_1);

    // grab the min and max for the r
    rMin = masHelBprDim_1[0];
    rMax = masHelBprDim_1[masHelBprDimMax_1[0] - 1];

    // masHelBtpDim
    masReadMeshDimensions(fileNames[1], "dim3", &masHelBtpDimMax_1[0]);
    masHelBtpDim_1 = (float *)malloc(sizeof(float) * (int)(masHelBtpDimMax_1[0]));
    masReadMesh(fileNames[1], "dim3", &masHelBtpDim_1);

    // masHelBttDim
    masReadMeshDimensions(fileNames[1], "dim2", &masHelBttDimMax_1[0]);
    masHelBttDim_1 = (float *)malloc(sizeof(float) * (int)(masHelBttDimMax_1[0]));
    masReadMesh(fileNames[1], "dim2", &masHelBttDim_1);

    // masHelBtrDim
    masReadMeshDimensions(fileNames[1], "dim1", &masHelBtrDimMax_1[0]);
    masHelBtrDim_1 = (float *)malloc(sizeof(float) * (int)(masHelBtrDimMax_1[0]));
    masReadMesh(fileNames[1], "dim1", &masHelBtrDim_1);

    // grab the min and max for the r
    rTemp = masHelBtrDim_1[0];
    if (rTemp > rMin)
      rMin = rTemp;

    rTemp = masHelBtrDim_1[masHelBtrDimMax_1[0] - 1];
    if (rTemp < rMax)
      rMax = rTemp;

    // masHelBrpDim
    masReadMeshDimensions(fileNames[2], "dim3", &masHelBrpDimMax_1[0]);
    masHelBrpDim_1 = (float *)malloc(sizeof(float) * (int)(masHelBrpDimMax_1[0]));
    masReadMesh(fileNames[2], "dim3", &masHelBrpDim_1);

    // masHelBrtDim
    masReadMeshDimensions(fileNames[2], "dim2", &masHelBrtDimMax_1[0]);
    masHelBrtDim_1 = (float *)malloc(sizeof(float) * (int)(masHelBrtDimMax_1[0]));
    masReadMesh(fileNames[2], "dim2", &masHelBrtDim_1);

    // masHelBrrDim
    masReadMeshDimensions(fileNames[2], "dim1", &masHelBrrDimMax_1[0]);
    masHelBrrDim_1 = (float *)malloc(sizeof(float) * (int)(masHelBrrDimMax_1[0]));
    masReadMesh(fileNames[2], "dim1", &masHelBrrDim_1);

    // grab the min and max for the r
    rTemp = masHelBrrDim_1[0];
    if (rTemp > rMin)
      rMin = rTemp;

    rTemp = masHelBrrDim_1[masHelBrrDimMax_1[0] - 1];
    if (rTemp < rMax)
      rMax = rTemp;

    // masHelVppDim
    masReadMeshDimensions(fileNames[3], "dim3", &masHelVppDimMax_1[0]);
    masHelVppDim_1 = (float *)malloc(sizeof(float) * (int)(masHelVppDimMax_1[0]));
    masReadMesh(fileNames[3], "dim3", &masHelVppDim_1);

    // masHelVptDim
    masReadMeshDimensions(fileNames[3], "dim2", &masHelVptDimMax_1[0]);
    masHelVptDim_1 = (float *)malloc(sizeof(float) * (int)(masHelVptDimMax_1[0]));
    masReadMesh(fileNames[3], "dim2", &masHelVptDim_1);

    // masHelVprDim
    masReadMeshDimensions(fileNames[3], "dim1", &masHelVprDimMax_1[0]);
    masHelVprDim_1 = (float *)malloc(sizeof(float) * (int)(masHelVprDimMax_1[0]));
    masReadMesh(fileNames[3], "dim1", &masHelVprDim_1);

    // grab the min and max for the r
    rTemp = masHelVprDim_1[0];
    if (rTemp > rMin)
      rMin = rTemp;

    rTemp = masHelVprDim_1[masHelVprDimMax_1[0] - 1];
    if (rTemp < rMax)
      rMax = rTemp;

    // masHelVtpDim
    masReadMeshDimensions(fileNames[4], "dim3", &masHelVtpDimMax_1[0]);
    masHelVtpDim_1 = (float *)malloc(sizeof(float) * (int)(masHelVtpDimMax_1[0]));
    masReadMesh(fileNames[4], "dim3", &masHelVtpDim_1);

    // masHelVttDim
    masReadMeshDimensions(fileNames[4], "dim2", &masHelVttDimMax_1[0]);
    masHelVttDim_1 = (float *)malloc(sizeof(float) * (int)(masHelVttDimMax_1[0]));
    masReadMesh(fileNames[4], "dim2", &masHelVttDim_1);

    // masHelVtrDim
    masReadMeshDimensions(fileNames[4], "dim1", &masHelVtrDimMax_1[0]);
    masHelVtrDim_1 = (float *)malloc(sizeof(float) * (int)(masHelVtrDimMax_1[0]));
    masReadMesh(fileNames[4], "dim1", &masHelVtrDim_1);

    // grab the min and max for the r
    rTemp = masHelVtrDim_1[0];
    if (rTemp > rMin)
      rMin = rTemp;

    rTemp = masHelVtrDim_1[masHelVtrDimMax_1[0] - 1];
    if (rTemp < rMax)
      rMax = rTemp;

    // masHelVrpDim
    masReadMeshDimensions(fileNames[5], "dim3", &masHelVrpDimMax_1[0]);
    masHelVrpDim_1 = (float *)malloc(sizeof(float) * (int)(masHelVrpDimMax_1[0]));
    masReadMesh(fileNames[5], "dim3", &masHelVrpDim_1);

    // masHelVrtDim
    masReadMeshDimensions(fileNames[5], "dim2", &masHelVrtDimMax_1[0]);
    masHelVrtDim_1 = (float *)malloc(sizeof(float) * (int)(masHelVrtDimMax_1[0]));
    masReadMesh(fileNames[5], "dim2", &masHelVrtDim_1);

    // masHelVrrDim
    masReadMeshDimensions(fileNames[5], "dim1", &masHelVrrDimMax_1[0]);
    masHelVrrDim_1 = (float *)malloc(sizeof(float) * (int)(masHelVrrDimMax_1[0]));
    masReadMesh(fileNames[5], "dim1", &masHelVrrDim_1);

    // grab the min and max for the r
    rTemp = masHelVrrDim_1[0];
    if (rTemp > rMin)
      rMin = rTemp;

    rTemp = masHelVrrDim_1[masHelVrrDimMax_1[0] - 1];
    if (rTemp < rMax)
      rMax = rTemp;

    // masHelDpDim
    masReadMeshDimensions(fileNames[6], "dim3", &masHelDpDimMax_1[0]);
    masHelDpDim_1 = (float *)malloc(sizeof(float) * (int)(masHelDpDimMax_1[0]));
    masReadMesh(fileNames[6], "dim3", &masHelDpDim_1);

    // masHelDtDim
    masReadMeshDimensions(fileNames[6], "dim2", &masHelDtDimMax_1[0]);
    masHelDtDim_1 = (float *)malloc(sizeof(float) * (int)(masHelDtDimMax_1[0]));
    masReadMesh(fileNames[6], "dim2", &masHelDtDim_1);

    // masHelDrDim
    masReadMeshDimensions(fileNames[6], "dim1", &masHelDrDimMax_1[0]);
    masHelDrDim_1 = (float *)malloc(sizeof(float) * (int)(masHelDrDimMax_1[0]));
    masReadMesh(fileNames[6], "dim1", &masHelDrDim_1);

    // grab the min and max for the r
    rTemp = masHelDrDim_1[0];
    if (rTemp > rMin)
      rMin = rTemp;

    rTemp = masHelDrDim_1[masHelDrDimMax_1[0] - 1];
    if (rTemp < rMax)
      rMax = rTemp;
  }
  else
  {
    printf("Unknown state \n");
    exit(0);
  }

  // set rScale, and masHelRadialMin/Max
  config.masHelRadialMin = rMin;
  config.masHelRadialMax = rMax;

  timer_mas_io = timer_mas_io + (MPI_Wtime() - timer_tmp);
} /*-------- END masHelReadFieldIndex()  -----------------------------*/
  /*------------------------------------------------------------------*/

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/
/*--*/ void                     /*--*/
                                /*--*/
allocateMPIWindows(char *state) /*--*/
                                /*--*/
/*--                                                                 --*/
/*--This function allocates the storage for mas variables.           --*/
/*---------------------------------------------------------------------*/
{ /*-------------------------------------------------------------------*/
  MPI_Aint N_0;
  MPI_Aint N_1;
  MPI_Aint size;
  int disp_unit;

  if (strcmp(state, "state0") == 0)
  {

    if (mpi_rank_shared == 0)
    {
      N_0 = (int)masBprDimMax_0[0] * (int)masBptDimMax_0[0] * (int)masBppDimMax_0[0];
      MPI_Win_allocate_shared(N_0 * sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masBp_0, &masBp_0_win);

      N_0 = (int)masBtrDimMax_0[0] * (int)masBttDimMax_0[0] * (int)masBtpDimMax_0[0];
      MPI_Win_allocate_shared(N_0 * sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masBt_0, &masBt_0_win);

      N_0 = (int)masBrrDimMax_0[0] * (int)masBrtDimMax_0[0] * (int)masBrpDimMax_0[0];
      MPI_Win_allocate_shared(N_0 * sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masBr_0, &masBr_0_win);

      N_0 = (int)masVprDimMax_0[0] * (int)masVptDimMax_0[0] * (int)masVppDimMax_0[0];
      MPI_Win_allocate_shared(N_0 * sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masVp_0, &masVp_0_win);

      N_0 = (int)masVtrDimMax_0[0] * (int)masVttDimMax_0[0] * (int)masVtpDimMax_0[0];
      MPI_Win_allocate_shared(N_0 * sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masVt_0, &masVt_0_win);

      N_0 = (int)masVrrDimMax_0[0] * (int)masVrtDimMax_0[0] * (int)masVrpDimMax_0[0];
      MPI_Win_allocate_shared(N_0 * sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masVr_0, &masVr_0_win);

      N_0 = (int)masDrDimMax_0[0] * (int)masDtDimMax_0[0] * (int)masDpDimMax_0[0];
      MPI_Win_allocate_shared(N_0 * sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masD_0, &masD_0_win);
    }

    else
    {
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masBp_0, &masBp_0_win);
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masBt_0, &masBt_0_win);
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masBr_0, &masBr_0_win);
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masVp_0, &masVp_0_win);
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masVt_0, &masVt_0_win);
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masVr_0, &masVr_0_win);
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masD_0, &masD_0_win);
      MPI_Win_shared_query(masBp_0_win, 0, &size, &disp_unit, &masBp_0);
      MPI_Win_shared_query(masBt_0_win, 0, &size, &disp_unit, &masBt_0);
      MPI_Win_shared_query(masBr_0_win, 0, &size, &disp_unit, &masBr_0);
      MPI_Win_shared_query(masVp_0_win, 0, &size, &disp_unit, &masVp_0);
      MPI_Win_shared_query(masVt_0_win, 0, &size, &disp_unit, &masVt_0);
      MPI_Win_shared_query(masVr_0_win, 0, &size, &disp_unit, &masVr_0);
      MPI_Win_shared_query(masD_0_win, 0, &size, &disp_unit, &masD_0);
    }
    MPI_Win_lock_all(MPI_MODE_NOCHECK, masBp_0_win);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, masBt_0_win);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, masBr_0_win);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, masVp_0_win);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, masVt_0_win);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, masVr_0_win);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, masD_0_win);
  }
  else if (strcmp(state, "state1") == 0)
  {

    if (mpi_rank_shared == 0)
    {
      N_1 = (int)masBprDimMax_1[0] * (int)masBptDimMax_1[0] * (int)masBppDimMax_1[0];
      MPI_Win_allocate_shared(N_1 * sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masBp_1, &masBp_1_win);

      N_1 = (int)masBtrDimMax_1[0] * (int)masBttDimMax_1[0] * (int)masBtpDimMax_1[0];
      MPI_Win_allocate_shared(N_1 * sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masBt_1, &masBt_1_win);

      N_1 = (int)masBrrDimMax_1[0] * (int)masBrtDimMax_1[0] * (int)masBrpDimMax_1[0];
      MPI_Win_allocate_shared(N_1 * sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masBr_1, &masBr_1_win);

      N_1 = (int)masVprDimMax_1[0] * (int)masVptDimMax_1[0] * (int)masVppDimMax_1[0];
      MPI_Win_allocate_shared(N_1 * sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masVp_1, &masVp_1_win);

      N_1 = (int)masVtrDimMax_1[0] * (int)masVttDimMax_1[0] * (int)masVtpDimMax_1[0];
      MPI_Win_allocate_shared(N_1 * sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masVt_1, &masVt_1_win);

      N_1 = (int)masVrrDimMax_1[0] * (int)masVrtDimMax_1[0] * (int)masVrpDimMax_1[0];
      MPI_Win_allocate_shared(N_1 * sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masVr_1, &masVr_1_win);

      N_1 = (int)masDrDimMax_1[0] * (int)masDtDimMax_1[0] * (int)masDpDimMax_1[0];
      MPI_Win_allocate_shared(N_1 * sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masD_1, &masD_1_win);
    }

    else
    {
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masBp_1, &masBp_1_win);
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masBt_1, &masBt_1_win);
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masBr_1, &masBr_1_win);
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masVp_1, &masVp_1_win);
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masVt_1, &masVt_1_win);
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masVr_1, &masVr_1_win);
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masD_1, &masD_1_win);
      MPI_Win_shared_query(masBp_1_win, 0, &size, &disp_unit, &masBp_1);
      MPI_Win_shared_query(masBt_1_win, 0, &size, &disp_unit, &masBt_1);
      MPI_Win_shared_query(masBr_1_win, 0, &size, &disp_unit, &masBr_1);
      MPI_Win_shared_query(masVp_1_win, 0, &size, &disp_unit, &masVp_1);
      MPI_Win_shared_query(masVt_1_win, 0, &size, &disp_unit, &masVt_1);
      MPI_Win_shared_query(masVr_1_win, 0, &size, &disp_unit, &masVr_1);
      MPI_Win_shared_query(masD_1_win, 0, &size, &disp_unit, &masD_1);
    }
    MPI_Win_lock_all(MPI_MODE_NOCHECK, masBp_1_win);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, masBt_1_win);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, masBr_1_win);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, masVp_1_win);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, masVt_1_win);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, masVr_1_win);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, masD_1_win);
  }
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/
/*--*/ void                        /*--*/
                                   /*--*/
allocateHelMPIWindows(char *state) /*--*/
                                   /*--*/
/*--                                                                 --*/
/*--This function allocates the storage for Hel MPIvariables.        --*/
/*---------------------------------------------------------------------*/
{ /*-------------------------------------------------------------------*/
  MPI_Aint N_0;
  MPI_Aint N_1;
  MPI_Aint size;
  int disp_unit;

  if (strcmp(state, "state0") == 0)
  {

    if (mpi_rank_shared == 0)
    {
      N_0 = (int)masHelBprDimMax_0[0] * (int)masHelBptDimMax_0[0] * (int)masHelBppDimMax_0[0];
      MPI_Win_allocate_shared(N_0 * sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masHelBp_0, &masHelBp_0_win);

      N_0 = (int)masHelBtrDimMax_0[0] * (int)masHelBttDimMax_0[0] * (int)masHelBtpDimMax_0[0];
      MPI_Win_allocate_shared(N_0 * sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masHelBt_0, &masHelBt_0_win);

      N_0 = (int)masHelBrrDimMax_0[0] * (int)masHelBrtDimMax_0[0] * (int)masHelBrpDimMax_0[0];
      MPI_Win_allocate_shared(N_0 * sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masHelBr_0, &masHelBr_0_win);

      N_0 = (int)masHelVprDimMax_0[0] * (int)masHelVptDimMax_0[0] * (int)masHelVppDimMax_0[0];
      MPI_Win_allocate_shared(N_0 * sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masHelVp_0, &masHelVp_0_win);

      N_0 = (int)masHelVtrDimMax_0[0] * (int)masHelVttDimMax_0[0] * (int)masHelVtpDimMax_0[0];
      MPI_Win_allocate_shared(N_0 * sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masHelVt_0, &masHelVt_0_win);

      N_0 = (int)masHelVrrDimMax_0[0] * (int)masHelVrtDimMax_0[0] * (int)masHelVrpDimMax_0[0];
      MPI_Win_allocate_shared(N_0 * sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masHelVr_0, &masHelVr_0_win);

      N_0 = (int)masHelDrDimMax_0[0] * (int)masHelDtDimMax_0[0] * (int)masHelDpDimMax_0[0];
      MPI_Win_allocate_shared(N_0 * sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masHelD_0, &masHelD_0_win);
    }

    else
    {
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masHelBp_0, &masHelBp_0_win);
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masHelBt_0, &masHelBt_0_win);
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masHelBr_0, &masHelBr_0_win);
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masHelVp_0, &masHelVp_0_win);
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masHelVt_0, &masHelVt_0_win);
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masHelVr_0, &masHelVr_0_win);
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masHelD_0, &masHelD_0_win);
      MPI_Win_shared_query(masHelBp_0_win, 0, &size, &disp_unit, &masHelBp_0);
      MPI_Win_shared_query(masHelBt_0_win, 0, &size, &disp_unit, &masHelBt_0);
      MPI_Win_shared_query(masHelBr_0_win, 0, &size, &disp_unit, &masHelBr_0);
      MPI_Win_shared_query(masHelVp_0_win, 0, &size, &disp_unit, &masHelVp_0);
      MPI_Win_shared_query(masHelVt_0_win, 0, &size, &disp_unit, &masHelVt_0);
      MPI_Win_shared_query(masHelVr_0_win, 0, &size, &disp_unit, &masHelVr_0);
      MPI_Win_shared_query(masHelD_0_win, 0, &size, &disp_unit, &masHelD_0);
    }

    MPI_Win_lock_all(MPI_MODE_NOCHECK, masHelBp_0_win);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, masHelBt_0_win);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, masHelBr_0_win);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, masHelVp_0_win);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, masHelVt_0_win);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, masHelVr_0_win);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, masHelD_0_win);
  }
  else if (strcmp(state, "state1") == 0)
  {

    if (mpi_rank_shared == 0)
    {

      N_1 = (int)masHelBprDimMax_1[0] * (int)masHelBptDimMax_1[0] * (int)masHelBppDimMax_1[0];
      MPI_Win_allocate_shared(N_1 * sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masHelBp_1, &masHelBp_1_win);

      N_1 = (int)masHelBtrDimMax_1[0] * (int)masHelBttDimMax_1[0] * (int)masHelBtpDimMax_1[0];
      MPI_Win_allocate_shared(N_1 * sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masHelBt_1, &masHelBt_1_win);

      N_1 = (int)masHelBrrDimMax_1[0] * (int)masHelBrtDimMax_1[0] * (int)masHelBrpDimMax_1[0];
      MPI_Win_allocate_shared(N_1 * sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masHelBr_1, &masHelBr_1_win);

      N_1 = (int)masHelVprDimMax_1[0] * (int)masHelVptDimMax_1[0] * (int)masHelVppDimMax_1[0];
      MPI_Win_allocate_shared(N_1 * sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masHelVp_1, &masHelVp_1_win);

      N_1 = (int)masHelVtrDimMax_1[0] * (int)masHelVttDimMax_1[0] * (int)masHelVtpDimMax_1[0];
      MPI_Win_allocate_shared(N_1 * sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masHelVt_1, &masHelVt_1_win);

      N_1 = (int)masHelVrrDimMax_1[0] * (int)masHelVrtDimMax_1[0] * (int)masHelVrpDimMax_1[0];
      MPI_Win_allocate_shared(N_1 * sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masHelVr_1, &masHelVr_1_win);

      N_1 = (int)masHelDrDimMax_1[0] * (int)masHelDtDimMax_1[0] * (int)masHelDpDimMax_1[0];
      MPI_Win_allocate_shared(N_1 * sizeof(float), sizeof(float), MPI_INFO_NULL, comm_shared, &masHelD_1, &masHelD_1_win);
    }

    else
    {

      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masHelBp_1, &masHelBp_1_win);
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masHelBt_1, &masHelBt_1_win);
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masHelBr_1, &masHelBr_1_win);
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masHelVp_1, &masHelVp_1_win);
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masHelVt_1, &masHelVt_1_win);
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masHelVr_1, &masHelVr_1_win);
      MPI_Win_allocate_shared(0, sizeof(float), MPI_INFO_NULL, comm_shared, &masHelD_1, &masHelD_1_win);
      MPI_Win_shared_query(masHelBp_1_win, 0, &size, &disp_unit, &masHelBp_1);
      MPI_Win_shared_query(masHelBt_1_win, 0, &size, &disp_unit, &masHelBt_1);
      MPI_Win_shared_query(masHelBr_1_win, 0, &size, &disp_unit, &masHelBr_1);
      MPI_Win_shared_query(masHelVp_1_win, 0, &size, &disp_unit, &masHelVp_1);
      MPI_Win_shared_query(masHelVt_1_win, 0, &size, &disp_unit, &masHelVt_1);
      MPI_Win_shared_query(masHelVr_1_win, 0, &size, &disp_unit, &masHelVr_1);
      MPI_Win_shared_query(masHelD_1_win, 0, &size, &disp_unit, &masHelD_1);
    }

    MPI_Win_lock_all(MPI_MODE_NOCHECK, masHelBp_1_win);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, masHelBt_1_win);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, masHelBr_1_win);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, masHelVp_1_win);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, masHelVt_1_win);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, masHelVr_1_win);
    MPI_Win_lock_all(MPI_MODE_NOCHECK, masHelD_1_win);
  }
}
/*--------------------------------------------------------------------*/

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/
/*--*/ int32_t                         /*--*/
                                     /*--*/
masReadBppDimTest(Index_t fileIndex) /*--*/
                                     /*--*/
/*--                                                                 --*/
/*--This function returns the dimension of Bpp mesh. Used to check   --*/
/*--any change in resolution                                          --*/
/*---------------------------------------------------------------------*/
{ /*-------------------------------------------------------------------*/

  char fileNames[7][MAX_STRING_SIZE];
  int32_t value_out[1];

  if (config.masDigits == 3)
  {
    sprintf(fileNames[0], "%sbp%03d%s", config.masDirectory, fileIndex + 1, file_extension);
  }
  else
  {
    sprintf(fileNames[0], "%sbp%06d%s", config.masDirectory, fileIndex + 1, file_extension);
  }

  // masBpp
  masReadMeshDimensions(fileNames[0], "dim3", &value_out[0]);

  return value_out[0];
}



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
#include <unistd.h>
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

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/
/*--*/ void                                                        /*--*/
                                                                   /*--*/
masReadMeshDimensions(char *fname, char *dsetname, int32_t *DimMax)/*--*/
/*--*/                                                             /*--*/
/*--                                                                 --*/
/*--Returns the number of elements in the 1D scale array "dsetname"  --*/
/*---------------------------------------------------------------------*/
{ /*-------------------------------------------------------------------*/

  hid_t file_id = -1;
  hid_t dataset_id = -1;
  hid_t dataspace_id = -1;
  hid_t file_dtype = -1;
  hsize_t size = -1;

  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT); // Open the file
  if (file_id < 0) {
    printf("ERROR: Cannot open file '%s'\n", fname);
    ERR(1);
  }

  dataset_id = H5Dopen2(file_id, dsetname, H5P_DEFAULT); // Open the dataset
  if (dataset_id < 0) {
    printf("ERROR: Cannot open dataset '%s' in '%s'\n", dsetname, fname);
    H5Fclose(file_id);
    ERR(1);
  }

  file_dtype = H5Dget_type(dataset_id);
  if (file_dtype < 0) {
    printf("ERROR: Could not get datatype for dataset '%s' in '%s'\n", dsetname, fname);
    H5Dclose(dataset_id);
    H5Fclose(file_id);
    ERR(1);
  }

  if (! H5Tequal(file_dtype, H5T_NATIVE_FLOAT)) {
    printf("ERROR: Dataset '%s' in '%s' is not H5T_NATIVE_FLOAT!\n", dsetname, fname);
    H5Tclose(file_dtype);
    H5Dclose(dataset_id);
    H5Fclose(file_id);
    ERR(1);
  }

  dataspace_id = H5Dget_space(dataset_id);
  if (dataspace_id < 0) {
    printf("ERROR: Cannot open dataspace '%s' in '%s'\n", dsetname, fname);
    H5Tclose(file_dtype);
    H5Dclose(dataset_id);
    H5Fclose(file_id);
    ERR(1);
  }

  size = H5Sget_simple_extent_npoints(dataspace_id);//, dim_sizes, NULL);
  if (size <= 0) {
    printf("ERROR: Size of mesh dimension for '%s' in '%s' is zero!\n", dsetname, fname);
    H5Sclose(dataspace_id);
    H5Tclose(file_dtype);
    H5Dclose(dataset_id);
    H5Fclose(file_id);
    ERR(1);
  }

  if (size > INT32_MAX) {
    printf("ERROR: Size of mesh dimension for '%s' in '%s' is too large for 32-bit integer: %lu elements\n", dsetname, fname, size);
    H5Sclose(dataspace_id);
    H5Tclose(file_dtype);
    H5Dclose(dataset_id);
    H5Fclose(file_id);
    ERR(1);
  }

  DimMax[0] = (int32_t)size;

  /* Clean up. */
  H5Sclose(dataspace_id);
  H5Tclose(file_dtype);
  H5Dclose(dataset_id);
  H5Fclose(file_id);
}
/*----------------- END masReadMeshDimensions() ------------------------*/
/*--------------------------------------------------------------------*/



/*---------------------------------------------------------------------*/
/*---------------------------------------------------------------------*/
/*---------------------------------------------------------------------*/
/*--*/ void                                                        /*--*/
                                                                   /*--*/
masReadMesh(char *fname, char *dsetname, float *Dim[])             /*--*/
/*--                                                                 --*/
/*--Returns the 1D scale array "dsetname".                            --*/
/*---------------------------------------------------------------------*/
{ /*-------------------------------------------------------------------*/
  hid_t file_id = -1;
  hid_t dataset_id = -1;
  hid_t file_dtype = -1;
  herr_t status = -1;

  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT); // Open the file
  if (file_id < 0) {
    printf("ERROR: Cannot open file '%s'\n", fname);
    ERR(1);
  }

  dataset_id = H5Dopen2(file_id, dsetname, H5P_DEFAULT); // Open the dataset
  if (dataset_id < 0) {
    printf("ERROR: Cannot open dataset '%s' in '%s'\n", dsetname, fname);
    H5Fclose(file_id);
    ERR(1);
  }

  file_dtype = H5Dget_type(dataset_id);
  if (file_dtype < 0) {
    printf("ERROR: Could not get datatype for dataset '%s' in '%s'\n", dsetname, fname);
    H5Dclose(dataset_id);
    H5Fclose(file_id);
    ERR(1);
  }

  if (! H5Tequal(file_dtype, H5T_NATIVE_FLOAT)) {
    printf("ERROR: Dataset '%s' in '%s' is not H5T_NATIVE_FLOAT!\n", dsetname, fname);
    H5Tclose(file_dtype);
    H5Dclose(dataset_id);
    H5Fclose(file_id);
    ERR(1);
  }

  status = H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, *Dim);
  if (status < 0) {
    printf("ERROR: Could not read dataset '%s' in '%s'\n", dsetname, fname);
    H5Tclose(file_dtype);
    H5Dclose(dataset_id);
    H5Fclose(file_id);
    ERR(1);
  }

  /* Clean up. */
  H5Tclose(file_dtype);
  H5Dclose(dataset_id);
  H5Fclose(file_id);
}
/*----------------- END masReadMesh() ----------------------------------*/
/*----------------------------------------------------------------------*/


/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
/*--*/ void                                                          /*--*/
                                                                     /*--*/
masReadDatafromFile(char *fname, float *buf[])                       /*--*/
/*--                                                                   --*/
/*--This function reads MAS 3D data stored in "Data"                   --*/
/*-----------------------------------------------------------------------*/
{ /*---------------------------------------------------------------------*/

  hid_t file_id = -1;
  hid_t dataset_id = -1;
  hid_t file_dtype = -1;
  herr_t status = -1;

  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT); // Open the file
  if (file_id < 0) {
    printf("ERROR: Cannot open file '%s'\n", fname);
    ERR(1);
  }

  dataset_id = H5Dopen2(file_id, "Data", H5P_DEFAULT); // Open the dataset
  if (dataset_id < 0) {
    printf("ERROR: Cannot open dataset Data in '%s'\n", fname);
    H5Fclose(file_id);
    ERR(1);
  }

  file_dtype = H5Dget_type(dataset_id);
  if (file_dtype < 0) {
    printf("ERROR: Could not get datatype for dataset Data in '%s'\n", fname);
    H5Dclose(dataset_id);
    H5Fclose(file_id);
    ERR(1);
  }

  if (! H5Tequal(file_dtype, H5T_NATIVE_FLOAT)) {
    printf("ERROR: Dataset Data in '%s' is not H5T_NATIVE_FLOAT!\n", fname);
    H5Tclose(file_dtype);
    H5Dclose(dataset_id);
    H5Fclose(file_id);
    ERR(1);
  }

  status = H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, *buf);
  if (status < 0) {
    printf("ERROR: Could not read Data for file '%s'\n", fname);
    H5Tclose(file_dtype);
    H5Dclose(dataset_id);
    H5Fclose(file_id);
    ERR(1);
  }

  /* Clean up. */
  H5Tclose(file_dtype);
  H5Dclose(dataset_id);
  H5Fclose(file_id);
}
/*----------------- END masReadDatafromfile() ----------------------------------*/
/*------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*--*/ void                                                         /*--*/
                                                                    /*--*/
masDatafile_type()                                                  /*--*/
                                                                    /*--*/
/*--                                                                  --*/
/*----------------------------------------------------------------------*/
{ /*--------------------------------------------------------------------*/

  char fileNames[7][MAX_STRING_SIZE];

  if (config.masDigits == 3)
  {
    sprintf(fileNames[0], "%sbp%03d", config.masDirectory, 1);
  }
  else
  {
    sprintf(fileNames[0], "%sbp%06d", config.masDirectory, 1);
  }

  if (access(strcat(fileNames[0],".h5"), F_OK) == 0) {
    printf("HDF5 mas datafiles detected \n");
    strncpy(file_extension,".h5",strlen(".h5")+1);
  }
  else{
    printf("No valid mas datafiles detected \n");
    exit(1);
  }
}

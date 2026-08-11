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

#include "global.h"
#include "configuration.h"
#include "masInterp.h"
#include "geometry.h"
#include "readMAS.h"
#include "mpiInit.h"
#include "error.h"
#include "simCore.h"
#include "flow.h"

masNode_t masNode;
Index_t unwindPhiOffset;

/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
/*--*/ void                                                          /*--*/
/*--*/ masGetNode(SphVec_t position, Node_t node)                    /*--*/
//     Interpolate MHD quantities to position.
//     Interpolating factors s_cor and s_hel computed in masGetInterpData()
//     to interpolate to desired time.  If this is called before
//     masGetInterpData() in main loop, the MHD quantities are
//     at time t_global (for example in RK4 in moveNodes).
//     If called after, the MHD quantities are at time t_global+dt.
//     The quantities are stored in the global "masNode" structure.
/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
{

  if (position.r <= (config.masRadialMax * RSAU))
  {

    masTriLinear(position,  masBp_0, masBt_0, masBr_0,
                            masVp_0, masVt_0, masVr_0,
                            masD_0,
                            masBp_1, masBt_1, masBr_1,
                            masVp_1, masVt_1, masVr_1,
                            masD_1,  s_cor);

  } else if ( (config.masHelCouple > 0) && (position.r <= (config.masHelRadialMax * RSAU)) ) {

    masHelTriLinear(position, masHelBp_0, masHelBt_0, masHelBr_0,
                              masHelVp_0, masHelVt_0, masHelVr_0,
                              masHelD_0,
                              masHelBp_1, masHelBt_1, masHelBr_1,
                              masHelVp_1, masHelVt_1, masHelVr_1,
                              masHelD_1,  s_hel);

  } else {

    masWind(node);

  }

}/*----------- END masGetNode() ------------------------------------*/
/*------------------------------------------------------------------*/


/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
/*--*/ int          /*---------------------------------------------------*/
/*--*/ masBinarySearch(float *A, float key, int imin, int imax)   /* ---*/
/*- does a binary search and returns the index                          -*/
/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
{
  int imid;

  while (imax != (imin + 1))
  {
    imid = (imin + imax) / 2;

    if (A[imid] < key) {

      imin = imid;

    } else {

      imax = imid;

    }
  }

  //if ((A[imin] > key) || (A[imax] < key)) {
  //  printf("Binary Search Fail!!!\n");
  //  printf("A[%i]:%f\tkry:%f\tA[%i]:%f\n", imin, A[imin], key, imax, A[imax]);
  //}

  if ((fabsf(A[imax] - key)) < (fabsf(A[imin] - key))) {

    return imax;

  } else {

    return imin;

  }

}
/*----------- END masBinarySearch() --------------------------------*/

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/
/*--*/ Scalar_t         /*------------------------------------------------------*/
/*--*/ masTriLinearBinarySearch(float *A, Scalar_t key, int *imin, int *imax,/**/
/*--*/                          int lower, int upper)                         /**/
/*- does a binary search and returns alpha, imin, and imax                     -*/
/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/
{
  int imid;
  Scalar_t returnValue;

  while (upper != (lower + 1))
  {
    imid = (lower + upper) / 2;

    if (A[imid] < key) {

      lower = imid;

    } else {

      upper = imid;

    }
  }

  *imin = lower;
  *imax = upper;

  //if ((A[*imin] > key) || (A[*imax] < key)) {
  //  printf("Binary Search Fail!!!\n");
  //  printf("A[%i]:%f\tkey:%f\tA[%i]:%f\n", lower, A[lower], key, upper, A[upper]);
  //}

  if (A[*imin] > key)
  {

    returnValue = 0.0;

  }
  else if (A[*imax] < key)
  {

    returnValue = 1.0;

  }
  else
  {

    returnValue = ( (key - A[lower]) / (A[upper] - A[lower]) );

  }

  return returnValue;

}
/*----------- END masTriLinearBinarySearch() -------------------------------*/


/*---------------------------------------------------------------------*/
/*---------------------------------------------------------------------*/
/*--*/ Scalar_t     /*-------------------------------------------------*/
/*--*/ masInterpolate(float V[], int r0, int r1, int t0, int t1,    /*-*/
/*--*/                int p0, int p1, Scalar_t rd, Scalar_t td,     /*-*/
/*--*/                Scalar_t pd, int rDimMax, int tDimMax)        /*-*/
/*--   Trilinearly interpolates and sends back the value              -*/
/*---------------------------------------------------------------------*/
/*---------------------------------------------------------------------*/
{
  #define Vindex(r,t,p) ((r) + rDimMax * (t) + rDimMax * tDimMax * (p))

  Scalar_t c[2][2];
  Scalar_t c0, c1;

  c[0][0] = V[Vindex(r0, t0, p0)] * (1.0 - rd) + V[Vindex(r1, t0, p0)] * rd;
  c[1][0] = V[Vindex(r0, t0, p1)] * (1.0 - rd) + V[Vindex(r1, t0, p1)] * rd;
  c[0][1] = V[Vindex(r0, t1, p0)] * (1.0 - rd) + V[Vindex(r1, t1, p0)] * rd;
  c[1][1] = V[Vindex(r0, t1, p1)] * (1.0 - rd) + V[Vindex(r1, t1, p1)] * rd;

  c0 = c[0][0] * (1.0 - pd) + c[1][0] * pd;
  c1 = c[0][1] * (1.0 - pd) + c[1][1] * pd;

  return (c0 * (1.0 - td) + c1 * td);

}
/*----------- END masInterpolate() ---------------------------------*/


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*--*/  SphVec_t                                                          /*--*/
/*--*/  fetchMasB( SphVec_t r,                                            /*--*/
/*--*/             float masBp0[], float masBt0[], float masBr0[],        /*--*/
/*--*/             float masBp1[], float masBt1[], float masBr1[],        /*--*/
/*--*/             Scalar_t s )                                           /*--*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
{

  int r0_0, r1_0, t0_0, t1_0, p0_0, p1_0;

  int r0_1, r1_1, t0_1, t1_1, p0_1, p1_1;

  Scalar_t rr_0, rd_0, td_0, pd_0;

  Scalar_t rd_1, td_1, pd_1;

  Scalar_t Bscale = 2.2068908 / MHD_B_NORM;

  SphVec_t B;

  if ((r.r < config.masRadialMax) && (r.r > config.masRadialMin)) {

    //Br
    rd_0 = masTriLinearBinarySearch(masBrrDim_0, r.r, &r0_0, &r1_0, masDimMin_0[0], masBrrDimMax_0[0]-1);
    td_0 = masTriLinearBinarySearch(masBrtDim_0, r.theta, &t0_0, &t1_0, masDimMin_0[0], masBrtDimMax_0[0]-1);
    pd_0 = masTriLinearBinarySearch(masBrpDim_0, r.phi, &p0_0, &p1_0, masDimMin_0[0], masBrpDimMax_0[0]-1);

    rd_1 = masTriLinearBinarySearch(masBrrDim_1, r.r, &r0_1, &r1_1, masDimMin_1[0], masBrrDimMax_1[0]-1);
    td_1 = masTriLinearBinarySearch(masBrtDim_1, r.theta, &t0_1, &t1_1, masDimMin_1[0], masBrtDimMax_1[0]-1);
    pd_1 = masTriLinearBinarySearch(masBrpDim_1, r.phi, &p0_1, &p1_1, masDimMin_1[0], masBrpDimMax_1[0]-1);

    B.r = ((1.0 - s) * masInterpolate(masBr0,
                                      r0_0, r1_0, t0_0, t1_0, p0_0, p1_0,
                                      rd_0, td_0, pd_0,
                                      masBrrDimMax_0[0], masBrtDimMax_0[0]) +
                   s * masInterpolate(masBr1,
                                      r0_1, r1_1, t0_1, t1_1, p0_1, p1_1,
                                      rd_1, td_1, pd_1,
                                      masBrrDimMax_1[0], masBrtDimMax_1[0])) * Bscale;

    //Bt
    rd_0 = masTriLinearBinarySearch(masBtrDim_0, r.r, &r0_0, &r1_0, masDimMin_0[0], masBtrDimMax_0[0]-1);
    td_0 = masTriLinearBinarySearch(masBttDim_0, r.theta, &t0_0, &t1_0, masDimMin_0[0], masBttDimMax_0[0]-1);
    pd_0 = masTriLinearBinarySearch(masBtpDim_0, r.phi, &p0_0, &p1_0, masDimMin_0[0], masBtpDimMax_0[0]-1);

    rd_1 = masTriLinearBinarySearch(masBtrDim_1, r.r, &r0_1, &r1_1, masDimMin_1[0], masBtrDimMax_1[0]-1);
    td_1 = masTriLinearBinarySearch(masBttDim_1, r.theta, &t0_1, &t1_1, masDimMin_1[0], masBttDimMax_1[0]-1);
    pd_1 = masTriLinearBinarySearch(masBtpDim_1, r.phi, &p0_1, &p1_1, masDimMin_1[0], masBtpDimMax_1[0]-1);

    B.theta = ((1.0 - s) * masInterpolate(masBt0,
                                          r0_0, r1_0, t0_0, t1_0, p0_0, p1_0,
                                          rd_0, td_0, pd_0,
                                          masBtrDimMax_0[0], masBttDimMax_0[0]) +
                       s * masInterpolate(masBt1,
                                          r0_1, r1_1, t0_1, t1_1, p0_1, p1_1,
                                          rd_1, td_1, pd_1,
                                          masBtrDimMax_1[0], masBttDimMax_1[0])) * Bscale;

    //Bp
    rd_0 = masTriLinearBinarySearch(masBprDim_0, r.r, &r0_0, &r1_0, masDimMin_0[0], masBprDimMax_0[0]-1);
    td_0 = masTriLinearBinarySearch(masBptDim_0, r.theta, &t0_0, &t1_0, masDimMin_0[0], masBptDimMax_0[0]-1);
    pd_0 = masTriLinearBinarySearch(masBppDim_0, r.phi, &p0_0, &p1_0, masDimMin_0[0], masBppDimMax_0[0]-1);

    rd_1 = masTriLinearBinarySearch(masBprDim_1, r.r, &r0_1, &r1_1, masDimMin_1[0], masBprDimMax_1[0]-1);
    td_1 = masTriLinearBinarySearch(masBptDim_1, r.theta, &t0_1, &t1_1, masDimMin_1[0], masBptDimMax_1[0]-1);
    pd_1 = masTriLinearBinarySearch(masBppDim_1, r.phi, &p0_1, &p1_1, masDimMin_1[0], masBppDimMax_1[0]-1);

    B.phi = ((1.0 - s) * masInterpolate(masBp0,
                                        r0_0, r1_0, t0_0, t1_0, p0_0, p1_0,
                                        rd_0, td_0, pd_0,
                                        masBprDimMax_0[0], masBptDimMax_0[0]) +
                     s * masInterpolate(masBp1,
                                        r0_1, r1_1, t0_1, t1_1, p0_1, p1_1,
                                        rd_1, td_1, pd_1,
                                        masBprDimMax_1[0], masBptDimMax_1[0])) * Bscale;

  } else {

    rr_0 = r.r * RSAU;

    B.r = config.mhdBsAu / (rr_0 * rr_0);

    B.theta = 0.0;

    B.phi = -1.0 * rr_0 * (B.r) * (config.omegaSun / (config.mhdUs + VERYSMALL) ) * sin(r.theta);

  }

  return B;

}
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*--*/ SphVec_t masCurlBoverB2( SphVec_t r,                               /*--*/
/*--*/                    float Bp0[], float Bt0[], float Br0[],          /*--*/
/*--*/                    float Bp1[], float Bt1[], float Br1[],          /*--*/
/*--*/                    int bp_r0_0, int bp_r1_0,                       /*--*/
/*--*/                    int bp_t0_0, int bp_t1_0,                       /*--*/
/*--*/                    int bp_p0_0, int bp_p1_0,                       /*--*/
/*--*/                    int bt_r0_0, int bt_r1_0,                       /*--*/
/*--*/                    int bt_t0_0, int bt_t1_0,                       /*--*/
/*--*/                    int bt_p0_0, int bt_p1_0,                       /*--*/
/*--*/                    int br_r0_0, int br_r1_0,                       /*--*/
/*--*/                    int br_t0_0, int br_t1_0,                       /*--*/
/*--*/                    int br_p0_0, int br_p1_0,                       /*--*/
/*--*/                    int bp_r0_1, int bp_r1_1,                       /*--*/
/*--*/                    int bp_t0_1, int bp_t1_1,                       /*--*/
/*--*/                    int bp_p0_1, int bp_p1_1,                       /*--*/
/*--*/                    int bt_r0_1, int bt_r1_1,                       /*--*/
/*--*/                    int bt_t0_1, int bt_t1_1,                       /*--*/
/*--*/                    int bt_p0_1, int bt_p1_1,                       /*--*/
/*--*/                    int br_r0_1, int br_r1_1,                       /*--*/
/*--*/                    int br_t0_1, int br_t1_1,                       /*--*/
/*--*/                    int br_p0_1, int br_p1_1,                       /*--*/
/*--*/                    Scalar_t s)                                     /*--*/
/*--*/                                                                    /*--*/
/*--*/                                                                    /*--*/
/*--     Calculate the curl of B/B^2 for use in the drift velocity          --*/
/*--     NOTE! This only works with Coronal MAS domain - Helio needs dev.   --*/
/*--                                                                        --*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
{

  Scalar_t phiCellWidth, phiCellWidthTemp;
  Scalar_t thetaCellWidth, thetaCellWidthTemp;
  Scalar_t rCellWidth, rCellWidthTemp;

  SphVec_t rMinus, rPlus, thetaMinus, thetaPlus, phiMinus, phiPlus;
  SphVec_t B_rm, B_rp, B_tm, B_tp, B_pm, B_pp;

  Scalar_t Bmag_rm, Bmag_rp, Bmag_tm, Bmag_tp, Bmag_pm, Bmag_pp;

  Scalar_t rAU, rCellWidthAU, A, B;

  SphVec_t curl;

  // Find the maximum cell width across both resolutions.
  // It is necessary to check each dimension of each component of the field
  // due to how the MAS field files are written out.

  // Bp - resolution 0
  rCellWidth = masBprDim_0[bp_r1_0] - masBprDim_0[bp_r0_0];
  thetaCellWidth = masBptDim_0[bp_t1_0] - masBptDim_0[bp_t0_0];
  phiCellWidth = masBppDim_0[bp_p1_0] - masBppDim_0[bp_p0_0];

  // Bt - resolution 0
  rCellWidthTemp = masBtrDim_0[bt_r1_0] - masBtrDim_0[bt_r0_0];
  thetaCellWidthTemp = masBttDim_0[bt_t1_0] - masBttDim_0[bt_t0_0];
  phiCellWidthTemp = masBtpDim_0[bt_p1_0] - masBtpDim_0[bt_p0_0];

  if (rCellWidthTemp > rCellWidth)
    rCellWidth = rCellWidthTemp;

  if (thetaCellWidthTemp > thetaCellWidth)
    thetaCellWidth = thetaCellWidthTemp;

  if (phiCellWidthTemp > phiCellWidth)
    phiCellWidth = phiCellWidthTemp;

  // Br - resolution 0
  rCellWidthTemp = masBrrDim_0[br_r1_0] - masBrrDim_0[br_r0_0];
  thetaCellWidthTemp = masBrtDim_0[br_t1_0] - masBrtDim_0[br_t0_0];
  phiCellWidthTemp = masBrpDim_0[br_p1_0] - masBrpDim_0[br_p0_0];

  if (rCellWidthTemp > rCellWidth)
    rCellWidth = rCellWidthTemp;

  if (thetaCellWidthTemp > thetaCellWidth)
    thetaCellWidth = thetaCellWidthTemp;

  if (phiCellWidthTemp > phiCellWidth)
    phiCellWidth = phiCellWidthTemp;

  // Bp - resolution 1
  rCellWidthTemp = masBprDim_1[bp_r1_1] - masBprDim_1[bp_r0_1];
  thetaCellWidthTemp = masBptDim_1[bp_t1_1] - masBptDim_1[bp_t0_1];
  phiCellWidthTemp = masBppDim_1[bp_p1_1] - masBppDim_1[bp_p0_1];

  if (rCellWidthTemp > rCellWidth)
    rCellWidth = rCellWidthTemp;

  if (thetaCellWidthTemp > thetaCellWidth)
    thetaCellWidth = thetaCellWidthTemp;

  if (phiCellWidthTemp > phiCellWidth)
    phiCellWidth = phiCellWidthTemp;

  // Bt - resolution 1
  rCellWidthTemp = masBtrDim_1[bt_r1_1] - masBtrDim_1[bt_r0_1];
  thetaCellWidthTemp = masBttDim_1[bt_t1_1] - masBttDim_1[bt_t0_1];
  phiCellWidthTemp = masBtpDim_1[bt_p1_1] - masBtpDim_1[bt_p0_1];

  if (rCellWidthTemp > rCellWidth)
    rCellWidth = rCellWidthTemp;

  if (thetaCellWidthTemp > thetaCellWidth)
    thetaCellWidth = thetaCellWidthTemp;

  if (phiCellWidthTemp > phiCellWidth)
    phiCellWidth = phiCellWidthTemp;

  // Br - resolution 1
  rCellWidthTemp = masBrrDim_1[br_r1_1] - masBrrDim_1[br_r0_1];
  thetaCellWidthTemp = masBrtDim_1[br_t1_1] - masBrtDim_1[br_t0_1];
  phiCellWidthTemp = masBrpDim_1[br_p1_1] - masBrpDim_1[br_p0_1];

  if (rCellWidthTemp > rCellWidth)
    rCellWidth = rCellWidthTemp;

  if (thetaCellWidthTemp > thetaCellWidth)
    thetaCellWidth = thetaCellWidthTemp;

  if (phiCellWidthTemp > phiCellWidth)
    phiCellWidth = phiCellWidthTemp;


  // Assign the spatial values for each delta direction and take into account
  // the boundaries for phi and theta
  rPlus.r = r.r + rCellWidth;
  rPlus.theta = r.theta;
  rPlus.phi = r.phi;

  rMinus.r = r.r - rCellWidth;
  rMinus.theta = r.theta;
  rMinus.phi = r.phi;

  thetaPlus.r = r.r;
  thetaPlus.theta = r.theta + thetaCellWidth;
  thetaPlus.phi = r.phi;

  if (thetaPlus.theta > PI)
  {
    thetaPlus.theta = 2.0 * PI - thetaPlus.theta;
    thetaPlus.phi += PI;

    if ( thetaPlus.phi > (2.0 * PI) )
      thetaPlus.phi -= (2.0 * PI);
  }

  thetaMinus.r = r.r;
  thetaMinus.theta = r.theta - thetaCellWidth;
  thetaMinus.phi = r.phi;

  if (thetaMinus.theta < 0.0)
  {
    thetaMinus.theta *= -1.0;
    thetaMinus.phi += PI;

    if ( thetaMinus.phi > (2.0 * PI) )
      thetaMinus.phi -= (2.0 * PI);
  }

  phiPlus.r = r.r;
  phiPlus.theta = r.theta;
  phiPlus.phi = r.phi + phiCellWidth;

  if ( phiPlus.phi > (2.0 * PI) )
    while ( phiPlus.phi > (2.0 * PI) ) phiPlus.phi = phiPlus.phi - (2.0 * PI);

  phiMinus.r = r.r;
  phiMinus.theta = r.theta;
  phiMinus.phi = r.phi - phiCellWidth;

  if ( phiMinus.phi < 0.0 )
    while ( phiMinus.phi < 0.0 ) phiMinus.phi = phiMinus.phi + (2.0 * PI);

  // Now that we have the spatial components, collect up the actual field values.
  B_rm = fetchMasB(rMinus,     Bp0, Bt0, Br0, Bp1, Bt1, Br1, s);
  B_rp = fetchMasB(rPlus,      Bp0, Bt0, Br0, Bp1, Bt1, Br1, s);

  B_tm = fetchMasB(thetaMinus, Bp0, Bt0, Br0, Bp1, Bt1, Br1, s);
  B_tp = fetchMasB(thetaPlus,  Bp0, Bt0, Br0, Bp1, Bt1, Br1, s);

  B_pm = fetchMasB(phiMinus,   Bp0, Bt0, Br0, Bp1, Bt1, Br1, s);
  B_pp = fetchMasB(phiPlus,    Bp0, Bt0, Br0, Bp1, Bt1, Br1, s);

  // calculate the magnitudes
  Bmag_rm = sqrt( B_rm.r * B_rm.r + B_rm.theta * B_rm.theta + B_rm.phi * B_rm.phi);
  Bmag_rp = sqrt( B_rp.r * B_rp.r + B_rp.theta * B_rp.theta + B_rp.phi * B_rp.phi);

  Bmag_tm = sqrt( B_tm.r * B_tm.r + B_tm.theta * B_tm.theta + B_tm.phi * B_tm.phi);
  Bmag_tp = sqrt( B_tp.r * B_tp.r + B_tp.theta * B_tp.theta + B_tp.phi * B_tp.phi);

  Bmag_pm = sqrt( B_pm.r * B_pm.r + B_pm.theta * B_pm.theta + B_pm.phi * B_pm.phi);
  Bmag_pp = sqrt( B_pp.r * B_pp.r + B_pp.theta * B_pp.theta + B_pp.phi * B_pp.phi);


  // finally!  calculating the curl itself
  rAU = r.r * RSAU;
  rCellWidthAU = rCellWidth * RSAU;

  // r
  A = (0.5 / thetaCellWidth) * ( B_tp.phi * sin(thetaPlus.theta) / (Bmag_tp * Bmag_tp) -
                                 B_tm.phi * sin(thetaMinus.theta) / (Bmag_tm * Bmag_tm) );

  B = (0.5 / phiCellWidth) * ( B_pp.theta / (Bmag_pp * Bmag_pp) - B_pm.theta / (Bmag_pm * Bmag_pm) );

  curl.r = (1.0 / (rAU * sin(r.theta))) * (A - B);

  // theta
  A = (1.0 / sin(r.theta)) * (0.5 / phiCellWidth)
      * ( B_pp.r / (Bmag_pp * Bmag_pp) - B_pm.r / (Bmag_pm * Bmag_pm) );

  B = (0.5 / rCellWidthAU) * ( (rPlus.r  * RSAU) * B_rp.phi / (Bmag_rp * Bmag_rp) -
                               (rMinus.r * RSAU) * B_rm.phi / (Bmag_rm * Bmag_rm) );

  curl.theta = (1.0 / rAU) * (A - B);

  // phi
  A = (0.5 / rCellWidthAU) * ( (rPlus.r  * RSAU) * B_rp.theta / (Bmag_rp * Bmag_rp) -
                               (rMinus.r * RSAU) * B_rm.theta / (Bmag_rm * Bmag_rm) );

  B = (0.5 / thetaCellWidth) * ( B_tp.r / (Bmag_tp * Bmag_tp)  - B_tm.r / (Bmag_tm * Bmag_tm) );

  curl.phi = (1.0 / rAU) * (A - B);


  return curl;

}
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*--*/ void masTriLinear( SphVec_t position,                              /*--*/
/*--*/                    float masBp0[], float masBt0[], float masBr0[], /*--*/
/*--*/                    float masVp0[], float masVt0[], float masVr0[], /*--*/
/*--*/                    float masD0[],                                  /*--*/
/*--*/                    float masBp1[], float masBt1[], float masBr1[], /*--*/
/*--*/                    float masVp1[], float masVt1[], float masVr1[], /*--*/
/*--*/                    float masD1[],                                  /*--*/
/*--*/                    Scalar_t s)                                     /*--*/
/*--*/                                                                    /*--*/
/*--*/                                                                    /*--*/
/*--  Get the data from the nearest mas node and interpolate                --*/
/*--  in time.                                                              --*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
{

  SphVec_t r;

  int r0_0, r1_0, t0_0, t1_0, p0_0, p1_0;
  int bp_r0_0, bp_r1_0, bp_t0_0, bp_t1_0, bp_p0_0, bp_p1_0;
  int bt_r0_0, bt_r1_0, bt_t0_0, bt_t1_0, bt_p0_0, bt_p1_0;
  int br_r0_0, br_r1_0, br_t0_0, br_t1_0, br_p0_0, br_p1_0;

  int r0_1, r1_1, t0_1, t1_1, p0_1, p1_1;
  int bp_r0_1, bp_r1_1, bp_t0_1, bp_t1_1, bp_p0_1, bp_p1_1;
  int bt_r0_1, bt_r1_1, bt_t0_1, bt_t1_1, bt_p0_1, bt_p1_1;
  int br_r0_1, br_r1_1, br_t0_1, br_t1_1, br_p0_1, br_p1_1;

  Scalar_t rd_0, td_0, pd_0;
  Scalar_t rd_1, td_1, pd_1;

  // convert from units of AU to Solar Radius and leave the
  //   angles alone
  r.r = position.r / RSAU;
  r.theta = position.theta;

  // do the angular offset in phi
  r.phi = position.phi - phiOffset;

  if (r.phi < 0.0)
    while (r.phi < 0.0) r.phi += (2.0 * PI);

  if ( r.phi > (2.0 * PI) )
    while ( r.phi > (2.0 * PI) ) r.phi -= (2.0 * PI);

  if ( (r.phi < 0.0) || (r.phi > 2.0 * PI) ) {
    if (mpi_rank == 0) printf("WARNING: r.phi (%f) out of bounds\n", r.phi);
  }

  //position
  masNode.r.r = position.r; // keep in units of AU
  masNode.r.theta = position.theta;
  masNode.r.phi = position.phi;


  //Bp
  rd_0 = masTriLinearBinarySearch(masBprDim_0, r.r, &bp_r0_0, &bp_r1_0, masDimMin_0[0], masBprDimMax_0[0]-1);
  td_0 = masTriLinearBinarySearch(masBptDim_0, r.theta, &bp_t0_0, &bp_t1_0, masDimMin_0[0], masBptDimMax_0[0]-1);
  pd_0 = masTriLinearBinarySearch(masBppDim_0, r.phi, &bp_p0_0, &bp_p1_0, masDimMin_0[0], masBppDimMax_0[0]-1);

  rd_1 = masTriLinearBinarySearch(masBprDim_1, r.r, &bp_r0_1, &bp_r1_1, masDimMin_1[0], masBprDimMax_1[0]-1);
  td_1 = masTriLinearBinarySearch(masBptDim_1, r.theta, &bp_t0_1, &bp_t1_1, masDimMin_1[0], masBptDimMax_1[0]-1);
  pd_1 = masTriLinearBinarySearch(masBppDim_1, r.phi, &bp_p0_1, &bp_p1_1, masDimMin_1[0], masBppDimMax_1[0]-1);

  masNode.mhdB.phi = ((1.0 - s) * masInterpolate(masBp0,
                                                bp_r0_0, bp_r1_0, bp_t0_0, bp_t1_0, bp_p0_0, bp_p1_0,
                                                rd_0, td_0, pd_0,
                                                masBprDimMax_0[0], masBptDimMax_0[0]) +
                              s * masInterpolate(masBp1,
                                                bp_r0_1, bp_r1_1, bp_t0_1, bp_t1_1, bp_p0_1, bp_p1_1,
                                                rd_1, td_1, pd_1,
                                                masBprDimMax_1[0], masBptDimMax_1[0])) * MAS_B_CONVERT;

  //Bt
  rd_0 = masTriLinearBinarySearch(masBtrDim_0, r.r, &bt_r0_0, &bt_r1_0, masDimMin_0[0], masBtrDimMax_0[0]-1);
  td_0 = masTriLinearBinarySearch(masBttDim_0, r.theta, &bt_t0_0, &bt_t1_0, masDimMin_0[0], masBttDimMax_0[0]-1);
  pd_0 = masTriLinearBinarySearch(masBtpDim_0, r.phi, &bt_p0_0, &bt_p1_0, masDimMin_0[0], masBtpDimMax_0[0]-1);

  rd_1 = masTriLinearBinarySearch(masBtrDim_1, r.r, &bt_r0_1, &bt_r1_1, masDimMin_1[0], masBtrDimMax_1[0]-1);
  td_1 = masTriLinearBinarySearch(masBttDim_1, r.theta, &bt_t0_1, &bt_t1_1, masDimMin_1[0], masBttDimMax_1[0]-1);
  pd_1 = masTriLinearBinarySearch(masBtpDim_1, r.phi, &bt_p0_1, &bt_p1_1, masDimMin_1[0], masBtpDimMax_1[0]-1);

  masNode.mhdB.theta = ((1.0 - s) * masInterpolate(masBt0,
                                                  bt_r0_0, bt_r1_0, bt_t0_0, bt_t1_0, bt_p0_0, bt_p1_0,
                                                  rd_0, td_0, pd_0,
                                                  masBtrDimMax_0[0], masBttDimMax_0[0]) +
                                s * masInterpolate(masBt1,
                                                  bt_r0_1, bt_r1_1, bt_t0_1, bt_t1_1, bt_p0_1, bt_p1_1,
                                                  rd_1, td_1, pd_1,
                                                  masBtrDimMax_1[0], masBttDimMax_1[0])) * MAS_B_CONVERT;

  //Br
  rd_0 = masTriLinearBinarySearch(masBrrDim_0, r.r, &br_r0_0, &br_r1_0, masDimMin_0[0], masBrrDimMax_0[0]-1);
  td_0 = masTriLinearBinarySearch(masBrtDim_0, r.theta, &br_t0_0, &br_t1_0, masDimMin_0[0], masBrtDimMax_0[0]-1);
  pd_0 = masTriLinearBinarySearch(masBrpDim_0, r.phi, &br_p0_0, &br_p1_0, masDimMin_0[0], masBrpDimMax_0[0]-1);

  rd_1 = masTriLinearBinarySearch(masBrrDim_1, r.r, &br_r0_1, &br_r1_1, masDimMin_1[0], masBrrDimMax_1[0]-1);
  td_1 = masTriLinearBinarySearch(masBrtDim_1, r.theta, &br_t0_1, &br_t1_1, masDimMin_1[0], masBrtDimMax_1[0]-1);
  pd_1 = masTriLinearBinarySearch(masBrpDim_1, r.phi, &br_p0_1, &br_p1_1, masDimMin_1[0], masBrpDimMax_1[0]-1);

  masNode.mhdB.r = ((1.0 - s) * masInterpolate(masBr0,
                                              br_r0_0, br_r1_0, br_t0_0, br_t1_0, br_p0_0, br_p1_0,
                                              rd_0, td_0, pd_0,
                                              masBrrDimMax_0[0], masBrtDimMax_0[0]) +
                            s * masInterpolate(masBr1,
                                              br_r0_1, br_r1_1, br_t0_1, br_t1_1, br_p0_1, br_p1_1,
                                              rd_1, td_1, pd_1,
                                              masBrrDimMax_1[0], masBrtDimMax_1[0])) * MAS_B_CONVERT;


  //Vp
  rd_0 = masTriLinearBinarySearch(masVprDim_0, r.r, &r0_0, &r1_0, masDimMin_0[0], masVprDimMax_0[0]-1);
  td_0 = masTriLinearBinarySearch(masVptDim_0, r.theta, &t0_0, &t1_0, masDimMin_0[0], masVptDimMax_0[0]-1);
  pd_0 = masTriLinearBinarySearch(masVppDim_0, r.phi, &p0_0, &p1_0, masDimMin_0[0], masVppDimMax_0[0]-1);

  rd_1 = masTriLinearBinarySearch(masVprDim_1, r.r, &r0_1, &r1_1, masDimMin_1[0], masVprDimMax_1[0]-1);
  td_1 = masTriLinearBinarySearch(masVptDim_1, r.theta, &t0_1, &t1_1, masDimMin_1[0], masVptDimMax_1[0]-1);
  pd_1 = masTriLinearBinarySearch(masVppDim_1, r.phi, &p0_1, &p1_1, masDimMin_1[0], masVppDimMax_1[0]-1);


  masNode.mhdV.phi = ((1.0 - s) * masInterpolate(masVp0,
                                                r0_0, r1_0, t0_0, t1_0, p0_0, p1_0,
                                                rd_0, td_0, pd_0,
                                                masVprDimMax_0[0], masVptDimMax_0[0]) +
                              s * masInterpolate(masVp1,
                                                r0_1, r1_1, t0_1, t1_1, p0_1, p1_1,
                                                rd_1, td_1, pd_1,
                                                masVprDimMax_1[0], masVptDimMax_1[0])) * MAS_V_CONVERT;

  //Vt
  rd_0 = masTriLinearBinarySearch(masVtrDim_0, r.r, &r0_0, &r1_0, masDimMin_0[0], masVtrDimMax_0[0]-1);
  td_0 = masTriLinearBinarySearch(masVttDim_0, r.theta, &t0_0, &t1_0, masDimMin_0[0], masVttDimMax_0[0]-1);
  pd_0 = masTriLinearBinarySearch(masVtpDim_0, r.phi, &p0_0, &p1_0, masDimMin_0[0], masVtpDimMax_0[0]-1);

  rd_1 = masTriLinearBinarySearch(masVtrDim_1, r.r, &r0_1, &r1_1, masDimMin_1[0], masVtrDimMax_1[0]-1);
  td_1 = masTriLinearBinarySearch(masVttDim_1, r.theta, &t0_1, &t1_1, masDimMin_1[0], masVttDimMax_1[0]-1);
  pd_1 = masTriLinearBinarySearch(masVtpDim_1, r.phi, &p0_1, &p1_1, masDimMin_1[0], masVtpDimMax_1[0]-1);

  masNode.mhdV.theta = ((1.0 - s) * masInterpolate(masVt0,
                                                  r0_0, r1_0, t0_0, t1_0, p0_0, p1_0,
                                                  rd_0, td_0, pd_0,
                                                  masVtrDimMax_0[0], masVttDimMax_0[0]) +
                                s * masInterpolate(masVt1,
                                                  r0_1, r1_1, t0_1, t1_1, p0_1, p1_1,
                                                  rd_1, td_1, pd_1,
                                                  masVtrDimMax_1[0], masVttDimMax_1[0])) * MAS_V_CONVERT;

  //Vr
  rd_0 = masTriLinearBinarySearch(masVrrDim_0, r.r, &r0_0, &r1_0, masDimMin_0[0], masVrrDimMax_0[0]-1);
  td_0 = masTriLinearBinarySearch(masVrtDim_0, r.theta, &t0_0, &t1_0, masDimMin_0[0], masVrtDimMax_0[0]-1);
  pd_0 = masTriLinearBinarySearch(masVrpDim_0, r.phi, &p0_0, &p1_0, masDimMin_0[0], masVrpDimMax_0[0]-1);

  rd_1 = masTriLinearBinarySearch(masVrrDim_1, r.r, &r0_1, &r1_1, masDimMin_1[0], masVrrDimMax_1[0]-1);
  td_1 = masTriLinearBinarySearch(masVrtDim_1, r.theta, &t0_1, &t1_1, masDimMin_1[0], masVrtDimMax_1[0]-1);
  pd_1 = masTriLinearBinarySearch(masVrpDim_1, r.phi, &p0_1, &p1_1, masDimMin_1[0], masVrpDimMax_1[0]-1);


  masNode.mhdV.r = ((1.0 - s) * masInterpolate(masVr0,
                                              r0_0, r1_0, t0_0, t1_0, p0_0, p1_0,
                                              rd_0, td_0, pd_0,
                                              masVrrDimMax_0[0], masVrtDimMax_0[0]) +
                            s * masInterpolate(masVr1,
                                              r0_1, r1_1, t0_1, t1_1, p0_1, p1_1,
                                              rd_1, td_1, pd_1,
                                              masVrrDimMax_1[0], masVrtDimMax_1[0])) * MAS_V_CONVERT;

  // check for underflows in Vr and set to min acceptable radial flow
  if ( masNode.mhdV.r < (config.masVmin / C) ) masNode.mhdV.r = (config.masVmin / C);

  //D
  rd_0 = masTriLinearBinarySearch(masDrDim_0, r.r, &r0_0, &r1_0, masDimMin_0[0], masDrDimMax_0[0]-1);
  td_0 = masTriLinearBinarySearch(masDtDim_0, r.theta, &t0_0, &t1_0, masDimMin_0[0], masDtDimMax_0[0]-1);
  pd_0 = masTriLinearBinarySearch(masDpDim_0, r.phi, &p0_0, &p1_0, masDimMin_0[0], masDpDimMax_0[0]-1);

  rd_1 = masTriLinearBinarySearch(masDrDim_1, r.r, &r0_1, &r1_1, masDimMin_1[0], masDrDimMax_1[0]-1);
  td_1 = masTriLinearBinarySearch(masDtDim_1, r.theta, &t0_1, &t1_1, masDimMin_1[0], masDtDimMax_1[0]-1);
  pd_1 = masTriLinearBinarySearch(masDpDim_1, r.phi, &p0_1, &p1_1, masDimMin_1[0], masDpDimMax_1[0]-1);

  masNode.mhdD = ((1.0 - s) * masInterpolate(masD0,
                                            r0_0, r1_0, t0_0, t1_0, p0_0, p1_0,
                                            rd_0, td_0, pd_0,
                                            masDrDimMax_0[0], masDtDimMax_0[0]) +
                          s * masInterpolate(masD1,
                                            r0_1, r1_1, t0_1, t1_1, p0_1, p1_1,
                                            rd_1, td_1, pd_1,
                                            masDrDimMax_1[0], masDtDimMax_1[0])) * MAS_RHO_CONVERT;


  // calculate the curl of B/B^2 if using shell drift
  if (config.useDrift > 0)
    masNode.curlBoverB2 = masCurlBoverB2(r,
                                         masBp0, masBt0, masBr0,
                                         masBp1, masBt1, masBr1,
                                         bp_r0_0, bp_r1_0, bp_t0_0, bp_t1_0, bp_p0_0, bp_p1_0,
                                         bt_r0_0, bt_r1_0, bt_t0_0, bt_t1_0, bt_p0_0, bt_p1_0,
                                         br_r0_0, br_r1_0, br_t0_0, br_t1_0, br_p0_0, br_p1_0,
                                         bp_r0_1, bp_r1_1, bp_t0_1, bp_t1_1, bp_p0_1, bp_p1_1,
                                         bt_r0_1, bt_r1_1, bt_t0_1, bt_t1_1, bt_p0_1, bt_p1_1,
                                         br_r0_1, br_r1_1, br_t0_1, br_t1_1, br_p0_1, br_p1_1,
                                         s);

}
/*----------- END masTriLinear() --------------------------------*/


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/
/*--*/ void masHelTriLinear( SphVec_t position,                             /*--*/
/*--*/                      float masBp0[], float masBt0[], float masBr0[], /*--*/
/*--*/                      float masVp0[], float masVt0[], float masVr0[], /*--*/
/*--*/                      float masD0[],                                  /*--*/
/*--*/                      float masBp1[], float masBt1[], float masBr1[], /*--*/
/*--*/                      float masVp1[], float masVt1[], float masVr1[], /*--*/
/*--*/                      float masD1[],                                  /*--*/
/*--*/                      Scalar_t s)                                     /*--*/
/*--*/                                                                      /*--*/
/*--*/                                                                      /*--*/
/*--  Get the data from the nearest mas node and interpolate                  --*/
/*--  in time.                                                                --*/
/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/
{

  SphVec_t r;

  int r0_0, r1_0, t0_0, t1_0, p0_0, p1_0;
  int bp_r0_0, bp_r1_0, bp_t0_0, bp_t1_0, bp_p0_0, bp_p1_0;
  int bt_r0_0, bt_r1_0, bt_t0_0, bt_t1_0, bt_p0_0, bt_p1_0;
  int br_r0_0, br_r1_0, br_t0_0, br_t1_0, br_p0_0, br_p1_0;

  int r0_1, r1_1, t0_1, t1_1, p0_1, p1_1;
  int bp_r0_1, bp_r1_1, bp_t0_1, bp_t1_1, bp_p0_1, bp_p1_1;
  int bt_r0_1, bt_r1_1, bt_t0_1, bt_t1_1, bt_p0_1, bt_p1_1;
  int br_r0_1, br_r1_1, br_t0_1, br_t1_1, br_p0_1, br_p1_1;

  Scalar_t rd_0, td_0, pd_0, phiHelSeedOffset;
  Scalar_t rd_1, td_1, pd_1;
  // If we are seeding nodes, the heliosphere is treated as
  // a co-rotating frame.  Thus, we need to shift back at
  // the corona's phiOffset.  Once the simulation starts,
  // the heliosphere is in the intertial frame so we no longer
  // need an additional offset.
  if (simStarted == 0){
    phiHelSeedOffset = phiOffset;
  } else {
    phiHelSeedOffset = 0.0;
  }

  // convert from units of AU to Solar Radius and leave the
  //   angles alone
  r.r = position.r / RSAU;
  r.theta = position.theta;

  // Set phi offset for getting data in heliosphere.
  // This is typically a fixed value, read from the metadata of the hel mhd data.
  // During seeding of nodes, an additional shift is applied to keep the heliosphere
  // aligned with the rotating coronal data.
  r.phi = position.phi - phiHelOffset - phiHelSeedOffset;

  if (r.phi < 0.0)
    while (r.phi < 0.0) r.phi += (2.0 * PI);

  if ( r.phi > (2.0 * PI) )
    while ( r.phi > (2.0 * PI) ) r.phi -= (2.0 * PI);

  //Set position:
  masNode.r.r = position.r; // keep in units of AU
  masNode.r.theta = position.theta;
  masNode.r.phi = position.phi;

  //Bp
  rd_0 = masTriLinearBinarySearch(masHelBprDim_0, r.r, &bp_r0_0, &bp_r1_0, masDimMin_0[0], masHelBprDimMax_0[0]-1);
  td_0 = masTriLinearBinarySearch(masHelBptDim_0, r.theta, &bp_t0_0, &bp_t1_0, masDimMin_0[0], masHelBptDimMax_0[0]-1);
  pd_0 = masTriLinearBinarySearch(masHelBppDim_0, r.phi, &bp_p0_0, &bp_p1_0, masDimMin_0[0], masHelBppDimMax_0[0]-1);

  rd_1 = masTriLinearBinarySearch(masHelBprDim_1, r.r, &bp_r0_1, &bp_r1_1, masDimMin_1[0], masHelBprDimMax_1[0]-1);
  td_1 = masTriLinearBinarySearch(masHelBptDim_1, r.theta, &bp_t0_1, &bp_t1_1, masDimMin_1[0], masHelBptDimMax_1[0]-1);
  pd_1 = masTriLinearBinarySearch(masHelBppDim_1, r.phi, &bp_p0_1, &bp_p1_1, masDimMin_1[0], masHelBppDimMax_1[0]-1);

  masNode.mhdB.phi = ((1.0 - s) * masInterpolate(masBp0,
                                                bp_r0_0, bp_r1_0, bp_t0_0, bp_t1_0, bp_p0_0, bp_p1_0,
                                                rd_0, td_0, pd_0,
                                                masHelBprDimMax_0[0], masHelBptDimMax_0[0]) +
                              s * masInterpolate(masBp1,
                                                bp_r0_1, bp_r1_1, bp_t0_1, bp_t1_1, bp_p0_1, bp_p1_1,
                                                rd_1, td_1, pd_1,
                                                masHelBprDimMax_1[0], masHelBptDimMax_1[0])) * MAS_B_CONVERT;

  //Bt
  rd_0 = masTriLinearBinarySearch(masHelBtrDim_0, r.r, &bt_r0_0, &bt_r1_0, masDimMin_0[0], masHelBtrDimMax_0[0]-1);
  td_0 = masTriLinearBinarySearch(masHelBttDim_0, r.theta, &bt_t0_0, &bt_t1_0, masDimMin_0[0], masHelBttDimMax_0[0]-1);
  pd_0 = masTriLinearBinarySearch(masHelBtpDim_0, r.phi, &bt_p0_0, &bt_p1_0, masDimMin_0[0], masHelBtpDimMax_0[0]-1);

  rd_1 = masTriLinearBinarySearch(masHelBtrDim_1, r.r, &bt_r0_1, &bt_r1_1, masDimMin_1[0], masHelBtrDimMax_1[0]-1);
  td_1 = masTriLinearBinarySearch(masHelBttDim_1, r.theta, &bt_t0_1, &bt_t1_1, masDimMin_1[0], masHelBttDimMax_1[0]-1);
  pd_1 = masTriLinearBinarySearch(masHelBtpDim_1, r.phi, &bt_p0_1, &bt_p1_1, masDimMin_1[0], masHelBtpDimMax_1[0]-1);
  masNode.mhdB.theta = ((1.0 - s) * masInterpolate(masBt0,
                                                  bt_r0_0, bt_r1_0, bt_t0_0, bt_t1_0, bt_p0_0, bt_p1_0,
                                                  rd_0, td_0, pd_0,
                                                  masHelBtrDimMax_0[0], masHelBttDimMax_0[0]) +
                                s * masInterpolate(masBt1,
                                                  bt_r0_1, bt_r1_1, bt_t0_1, bt_t1_1, bt_p0_1, bt_p1_1,
                                                  rd_1, td_1, pd_1,
                                                  masHelBtrDimMax_1[0], masHelBttDimMax_1[0])) * MAS_B_CONVERT;

  //Br
  rd_0 = masTriLinearBinarySearch(masHelBrrDim_0, r.r, &br_r0_0, &br_r1_0, masDimMin_0[0], masHelBrrDimMax_0[0]-1);
  td_0 = masTriLinearBinarySearch(masHelBrtDim_0, r.theta, &br_t0_0, &br_t1_0, masDimMin_0[0], masHelBrtDimMax_0[0]-1);
  pd_0 = masTriLinearBinarySearch(masHelBrpDim_0, r.phi, &br_p0_0, &br_p1_0, masDimMin_0[0], masHelBrpDimMax_0[0]-1);

  rd_1 = masTriLinearBinarySearch(masHelBrrDim_1, r.r, &br_r0_1, &br_r1_1, masDimMin_1[0], masHelBrrDimMax_1[0]-1);
  td_1 = masTriLinearBinarySearch(masHelBrtDim_1, r.theta, &br_t0_1, &br_t1_1, masDimMin_1[0], masHelBrtDimMax_1[0]-1);
  pd_1 = masTriLinearBinarySearch(masHelBrpDim_1, r.phi, &br_p0_1, &br_p1_1, masDimMin_1[0], masHelBrpDimMax_1[0]-1);

  masNode.mhdB.r = ((1.0 - s) * masInterpolate(masBr0,
                                              br_r0_0, br_r1_0, br_t0_0, br_t1_0, br_p0_0, br_p1_0,
                                              rd_0, td_0, pd_0,
                                              masHelBrrDimMax_0[0], masHelBrtDimMax_0[0]) +
                            s * masInterpolate(masBr1,
                                              br_r0_1, br_r1_1, br_t0_1, br_t1_1, br_p0_1, br_p1_1,
                                              rd_1, td_1, pd_1,
                                              masHelBrrDimMax_1[0], masHelBrtDimMax_1[0])) * MAS_B_CONVERT;


  //Vp
  rd_0 = masTriLinearBinarySearch(masHelVprDim_0, r.r, &r0_0, &r1_0, masDimMin_0[0], masHelVprDimMax_0[0]-1);
  td_0 = masTriLinearBinarySearch(masHelVptDim_0, r.theta, &t0_0, &t1_0, masDimMin_0[0], masHelVptDimMax_0[0]-1);
  pd_0 = masTriLinearBinarySearch(masHelVppDim_0, r.phi, &p0_0, &p1_0, masDimMin_0[0], masHelVppDimMax_0[0]-1);

  rd_1 = masTriLinearBinarySearch(masHelVprDim_1, r.r, &r0_1, &r1_1, masDimMin_1[0], masHelVprDimMax_1[0]-1);
  td_1 = masTriLinearBinarySearch(masHelVptDim_1, r.theta, &t0_1, &t1_1, masDimMin_1[0], masHelVptDimMax_1[0]-1);
  pd_1 = masTriLinearBinarySearch(masHelVppDim_1, r.phi, &p0_1, &p1_1, masDimMin_1[0], masHelVppDimMax_1[0]-1);

  masNode.mhdV.phi = ((1.0 - s) * masInterpolate(masVp0,
                                                r0_0, r1_0, t0_0, t1_0, p0_0, p1_0,
                                                rd_0, td_0, pd_0,
                                                masHelVprDimMax_0[0], masHelVptDimMax_0[0]) +
                              s * masInterpolate(masVp1,
                                                r0_1, r1_1, t0_1, t1_1, p0_1, p1_1,
                                                rd_1, td_1, pd_1,
                                                masHelVprDimMax_1[0], masHelVptDimMax_1[0])) * MAS_V_CONVERT;

  //Vt
  rd_0 = masTriLinearBinarySearch(masHelVtrDim_0, r.r, &r0_0, &r1_0, masDimMin_0[0], masHelVtrDimMax_0[0]-1);
  td_0 = masTriLinearBinarySearch(masHelVttDim_0, r.theta, &t0_0, &t1_0, masDimMin_0[0], masHelVttDimMax_0[0]-1);
  pd_0 = masTriLinearBinarySearch(masHelVtpDim_0, r.phi, &p0_0, &p1_0, masDimMin_0[0], masHelVtpDimMax_0[0]-1);

  rd_1 = masTriLinearBinarySearch(masHelVtrDim_1, r.r, &r0_1, &r1_1, masDimMin_1[0], masHelVtrDimMax_1[0]-1);
  td_1 = masTriLinearBinarySearch(masHelVttDim_1, r.theta, &t0_1, &t1_1, masDimMin_1[0], masHelVttDimMax_1[0]-1);
  pd_1 = masTriLinearBinarySearch(masHelVtpDim_1, r.phi, &p0_1, &p1_1, masDimMin_1[0], masHelVtpDimMax_1[0]-1);

  masNode.mhdV.theta = ((1.0 - s) * masInterpolate(masVt0,
                                                  r0_0, r1_0, t0_0, t1_0, p0_0, p1_0,
                                                  rd_0, td_0, pd_0,
                                                  masHelVtrDimMax_0[0], masHelVttDimMax_0[0]) +
                                s * masInterpolate(masVt1,
                                                  r0_1, r1_1, t0_1, t1_1, p0_1, p1_1,
                                                  rd_1, td_1, pd_1,
                                                  masHelVtrDimMax_1[0], masHelVttDimMax_1[0])) * MAS_V_CONVERT;

  //Vr
  rd_0 = masTriLinearBinarySearch(masHelVrrDim_0, r.r, &r0_0, &r1_0, masDimMin_0[0], masHelVrrDimMax_0[0]-1);
  td_0 = masTriLinearBinarySearch(masHelVrtDim_0, r.theta, &t0_0, &t1_0, masDimMin_0[0], masHelVrtDimMax_0[0]-1);
  pd_0 = masTriLinearBinarySearch(masHelVrpDim_0, r.phi, &p0_0, &p1_0, masDimMin_0[0], masHelVrpDimMax_0[0]-1);

  rd_1 = masTriLinearBinarySearch(masHelVrrDim_1, r.r, &r0_1, &r1_1, masDimMin_1[0], masHelVrrDimMax_1[0]-1);
  td_1 = masTriLinearBinarySearch(masHelVrtDim_1, r.theta, &t0_1, &t1_1, masDimMin_1[0], masHelVrtDimMax_1[0]-1);
  pd_1 = masTriLinearBinarySearch(masHelVrpDim_1, r.phi, &p0_1, &p1_1, masDimMin_1[0], masHelVrpDimMax_1[0]-1);

  masNode.mhdV.r = ((1.0 - s) * masInterpolate(masVr0,
                                              r0_0, r1_0, t0_0, t1_0, p0_0, p1_0,
                                              rd_0, td_0, pd_0,
                                              masHelVrrDimMax_0[0], masHelVrtDimMax_0[0]) +
                            s * masInterpolate(masVr1,
                                              r0_1, r1_1, t0_1, t1_1, p0_1, p1_1,
                                              rd_1, td_1, pd_1,
                                              masHelVrrDimMax_1[0], masHelVrtDimMax_1[0])) * MAS_V_CONVERT;

  // check for underflows in Vr and set to min acceptable radial flow
  if ( masNode.mhdV.r < (config.masVmin / C) ) masNode.mhdV.r = (config.masVmin / C);

  //D
  rd_0 = masTriLinearBinarySearch(masHelDrDim_0, r.r, &r0_0, &r1_0, masDimMin_0[0], masHelDrDimMax_0[0]-1);
  td_0 = masTriLinearBinarySearch(masHelDtDim_0, r.theta, &t0_0, &t1_0, masDimMin_0[0], masHelDtDimMax_0[0]-1);
  pd_0 = masTriLinearBinarySearch(masHelDpDim_0, r.phi, &p0_0, &p1_0, masDimMin_0[0], masHelDpDimMax_0[0]-1);

  rd_1 = masTriLinearBinarySearch(masHelDrDim_1, r.r, &r0_1, &r1_1, masDimMin_1[0], masHelDrDimMax_1[0]-1);
  td_1 = masTriLinearBinarySearch(masHelDtDim_1, r.theta, &t0_1, &t1_1, masDimMin_1[0], masHelDtDimMax_1[0]-1);
  pd_1 = masTriLinearBinarySearch(masHelDpDim_1, r.phi, &p0_1, &p1_1, masDimMin_1[0], masHelDpDimMax_1[0]-1);

  masNode.mhdD = ((1.0 - s) * masInterpolate(masD0,
                                            r0_0, r1_0, t0_0, t1_0, p0_0, p1_0,
                                            rd_0, td_0, pd_0,
                                            masHelDrDimMax_0[0], masHelDtDimMax_0[0]) +
                          s * masInterpolate(masD1,
                                            r0_1, r1_1, t0_1, t1_1, p0_1, p1_1,
                                            rd_1, td_1, pd_1,
                                            masHelDrDimMax_1[0], masHelDtDimMax_1[0])) * MAS_RHO_CONVERT;

  // NOTE!  This needs to be implemented before using this on helio runs!

  // calculate the curl of B/B^2 if using shell drift
//  if (config.useDrift > 0)
//    masNode.curlBoverB2 = masHelCurlBoverB2(r,
//                                         masBp0, masBt0, masBr0,
//                                         masBp1, masBt1, masBr1,
//                                         bp_r0, bp_r1, bp_t0, bp_t1, bp_p0, bp_p1,
//                                         bt_r0, bt_r1, bt_t0, bt_t1, bt_p0, bt_p1,
//                                         br_r0, br_r1, br_t0, br_t1, br_p0, br_p1,
//                                         s);

}
/*----------- END masHelTriLinear() --------------------------------*/


/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/*--*/ SphVec_t fieldAlignedFlow(SphVec_t B)                    /*--*/
/*--*/                                                          /*--*/
/*--   Field aligned flow                                         --*/
/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
{

  SphVec_t V;
  Scalar_t Bmag = 0.0;

  if ( (fabs(B.r) * MHD_B_NORM) < .1 )
  {

    V.r     = config.parallelFlow / C;
    V.theta = 0.0;
    V.phi   = 0.0;

  }
  else if (B.r < 0.0)
  {

    B.r     *= -1.0;
    B.theta *= -1.0;
    B.phi   *= -1.0;

    Bmag = sqrt( B.r * B.r + B.theta * B.theta + B.phi * B.phi );

    V.r     = (config.parallelFlow / C) * B.r / Bmag;
    V.theta = (config.parallelFlow / C) * B.theta / Bmag;
    V.phi   = (config.parallelFlow / C) * B.phi / Bmag;

  }
  else
  {

    Bmag = sqrt( B.r * B.r + B.theta * B.theta + B.phi * B.phi );

    V.r     = (config.parallelFlow / C) * B.r / Bmag;
    V.theta = (config.parallelFlow / C) * B.theta / Bmag;
    V.phi   = (config.parallelFlow / C) * B.phi / Bmag;

  }

  return V;

}
/*----------- END fieldAlignedFlow() --------------------------------*/


/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/*--*/ void masWind(Node_t node)                                /*--*/
/*--*/                                                          /*--*/
/*--  Revert to the Parker wind model.                            --*/
/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
{
  Vec_t rOld;

  Scalar_t rmag, rOldmag, oneOverR, theta, thetaOld, br, bt, bp, vr ,vt ,vp, rho;

  rOld  = node.rOld;
  br    = node.mhdBr;
  bt    = node.mhdBtheta;
  bp    = node.mhdBphi;
  vr    = node.mhdVr;
  vt    = node.mhdVtheta;
  vp    = node.mhdVphi;
  rho   = node.mhdDensity;

  rOldmag = sqrt(rOld.x * rOld.x + rOld.y * rOld.y + rOld.z * rOld.z);
  rmag = node.rmag;

  oneOverR = rOldmag / rmag;

  // magnetic field
  masNode.mhdB.r     =  br * oneOverR * oneOverR;
  masNode.mhdB.theta =  bt * oneOverR;
  masNode.mhdB.phi   =  bp * oneOverR;

  // if we're co-rotating the inner solution then tack on the necessary phi component once outside.
  // If performing a heliospheric coupled run, this is not needed.
  // [RMC] Is this only needed for 'fake' corotation?
  // It seems we should NOT do this for true corotation...
  if ((config.masCorRotateFake > 0) && (config.masRotateSolution > 0) && (config.masHelCouple == 0))
  {
    theta = acos(node.r.z / rmag);
    thetaOld = acos(rOld.z / rOldmag);

    masNode.mhdB.phi +=
    ( (config.rScale * config.omegaSun / vr) * (rOldmag * br * sin(thetaOld) - rmag * masNode.mhdB.r * sin(theta)) );
  }

  // There is a slight issue when moving across the boundary (hitting a constant velocity) which is causing some
  // jagged bits along the node histories.  It isn't causing any problems but it
  // visually is annoying.  Fix when there is nothing critical going on.  Haha. - MG

  // velocity field
  masNode.mhdV.r     = vr;
  masNode.mhdV.theta = vt * oneOverR;
  masNode.mhdV.phi   = vp * oneOverR;

  // density
  masNode.mhdD = rho * oneOverR * oneOverR;

  // curlBoverB2
  if (config.useDrift > 0)
    masNode.curlBoverB2 = curlBoverB2(node.r, node.mhdVr, 0);

}
/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/

/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/*--*/ void rotateCoupledDomain( Scalar_t dt )                  /*--*/
/*--*/                                                          /*--*/
/*--  If inside the coupled (corotating) domain, rotate the nodes.--*/
/*--  This is equivalent to adding the proper +v_phi velocity     --*/
/*--  correction to the MHD.  This is (probably) more accurate.   --*/
/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
{

  Index_t face, row, col, shell, pObsIndex;

  Scalar_t rmag, dtPhiOffset;

  SphVec_t rSpherical, rOldSpherical;

  // Set the phi offset for this timestep
  dtPhiOffset = dt * config.omegaSun;

  // Set the phi offset for the MAS solution rotation
  phiOffset += dtPhiOffset;

  for (face = 0; face < NUM_FACES; face++)
  {
    for (row = 0; row < FACE_ROWS; row++)
    {
      for (col = 0; col < FACE_COLS; col++)
      {
        for (shell = INNER_ACTIVE_SHELL; shell < LOCAL_NUM_SHELLS; shell++)
        {

          // rmag in units of solar radius
          rmag = grid[idx_frcs(face,row,col,shell)].rmag * config.rScale / RSAU;

          if ( (rmag >= config.masRadialMin) && (rmag <= config.masRadialMax) )
          {

            // Apply the phi offset
            rSpherical = cartToSphPos(grid[idx_frcs(face,row,col,shell)].r);
            rSpherical.phi += dtPhiOffset;

            if (rSpherical.phi < 0.0)
              while (rSpherical.phi < 0.0) rSpherical.phi += (2.0 * PI);

            if ( rSpherical.phi > (2.0 * PI) )
              while ( rSpherical.phi > (2.0 * PI) ) rSpherical.phi -= (2.0 * PI);

            if ( (rSpherical.phi < 0.0) || (rSpherical.phi > 2.0 * PI) ) {
              if (mpi_rank == 0) printf("WARNING: rSpherical.phi (%f) out of bounds\n", rSpherical.phi);
            }

            rOldSpherical = cartToSphPos(grid[idx_frcs(face,row,col,shell)].rOld);
            rOldSpherical.phi += dtPhiOffset;

            if (rOldSpherical.phi < 0.0)
              while (rOldSpherical.phi < 0.0) rOldSpherical.phi += (2.0 * PI);

            if ( rOldSpherical.phi > (2.0 * PI) )
              while ( rOldSpherical.phi > (2.0 * PI) ) rOldSpherical.phi -= (2.0 * PI);

            if ( (rOldSpherical.phi < 0.0) || (rOldSpherical.phi > 2.0 * PI) ) {
              if (mpi_rank == 0) printf("WARNING: rOldSpherical.phi (%f) out of bounds\n", rOldSpherical.phi);
            }

            grid[idx_frcs(face,row,col,shell)].r = sphToCartPos(rSpherical);
            grid[idx_frcs(face,row,col,shell)].rOld = sphToCartPos(rOldSpherical);

          }
        }

      }

    }

  }

  // check for point observers living inside the coupled domain
  // RMC: this is probably broken for helio runs?
  for (pObsIndex = 0; pObsIndex < config.numObservers; pObsIndex++)
  {

    rmag = config.obsR[pObsIndex] / RSAU;

    // are we in the coupled domain?
    // if so, adjust the position through the entire run
    if ( (rmag > config.masRadialMin) && (rmag <= config.masRadialMax) )
    {

      config.obsPhi[pObsIndex] += dtPhiOffset;

      if (config.obsPhi[pObsIndex] < 0.0)
        while (config.obsPhi[pObsIndex] < 0.0) config.obsPhi[pObsIndex] += (2.0 * PI);

      if ( config.obsPhi[pObsIndex] > (2.0 * PI) )
        while ( config.obsPhi[pObsIndex] > (2.0 * PI) ) config.obsPhi[pObsIndex] -= (2.0 * PI);

      if ( (config.obsPhi[pObsIndex] < 0.0) || (config.obsPhi[pObsIndex] > 2.0 * PI) ) {
        if (mpi_rank == 0) printf("WARNING: config.obsPhi[pObsIndex] (%f) out of bounds\n", config.obsPhi[pObsIndex]);
      }

    }

    // are we outside of the coupled domain?
    // if so, adjust the position only during the equilibrium phase.
    if ( (rmag > config.masRadialMax) && (simStarted == 0) )
    {

      config.obsPhi[pObsIndex] += dtPhiOffset;

      if (config.obsPhi[pObsIndex] < 0.0)
        while (config.obsPhi[pObsIndex] < 0.0) config.obsPhi[pObsIndex] += (2.0 * PI);

      if ( config.obsPhi[pObsIndex] > (2.0 * PI) )
        while ( config.obsPhi[pObsIndex] > (2.0 * PI) ) config.obsPhi[pObsIndex] -= (2.0 * PI);

      if ( (config.obsPhi[pObsIndex] < 0.0) || (config.obsPhi[pObsIndex] > 2.0 * PI) ) {
        if (mpi_rank == 0) printf("WARNING: config.obsPhi[pObsIndex] (%f) out of bounds\n", config.obsPhi[pObsIndex]);
      }

    }


  }


}
/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/*--*/ void resetDomainOffset( void )                           /*--*/
/*--*/                                                          /*--*/
/*--  After the initial propogation of nodes, reset the offsets   --*/
/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
{

  Index_t face, row, col, shell, pObsIndex;
  SphVec_t rSpherical, rOldSpherical;

  for (face = 0; face < NUM_FACES; face++) {
    for (row = 0; row < FACE_ROWS; row++) {
      for (col = 0; col < FACE_COLS; col++) {
        for (shell = INNER_ACTIVE_SHELL; shell < LOCAL_NUM_SHELLS; shell++) {

          // reset the phi offset for both r and rOld
          rSpherical = cartToSphPos(grid[idx_frcs(face,row,col,shell)].r);
          rSpherical.phi -= phiOffset;

          if (rSpherical.phi < 0.0)
            while (rSpherical.phi < 0.0) rSpherical.phi += (2.0 * PI);

          if (rSpherical.phi > (2.0 * PI) )
            while (rSpherical.phi > (2.0 * PI) ) rSpherical.phi -= (2.0 * PI);

          if ( (rSpherical.phi < 0.0) || (rSpherical.phi > 2.0 * PI) ) {
            if (mpi_rank == 0) printf("WARNING: rSpherical.phi (%f) out of bounds\n", rSpherical.phi);
          }

          rOldSpherical = cartToSphPos(grid[idx_frcs(face,row,col,shell)].rOld);
          rOldSpherical.phi -= phiOffset;

          if (rOldSpherical.phi < 0.0)
            while (rOldSpherical.phi < 0.0) rOldSpherical.phi += (2.0 * PI);

          if (rOldSpherical.phi > (2.0 * PI) )
            while (rOldSpherical.phi > (2.0 * PI) ) rOldSpherical.phi -= (2.0 * PI);

          if ( (rOldSpherical.phi < 0.0) || (rOldSpherical.phi > 2.0 * PI) ) {
            if (mpi_rank == 0) printf("WARNING: rOldSpherical.phi (%f) out of bounds\n", rOldSpherical.phi);
          }

          grid[idx_frcs(face,row,col,shell)].r = sphToCartPos(rSpherical);
          grid[idx_frcs(face,row,col,shell)].rOld = sphToCartPos(rOldSpherical);

        }
      }
    }
  }

  // adjust the point observers
  for (pObsIndex = 0; pObsIndex < config.numObservers; pObsIndex++)
  {

      config.obsPhi[pObsIndex] -= phiOffset;

      if (config.obsPhi[pObsIndex] < 0.0)
        config.obsPhi[pObsIndex] += (2.0 * PI);

      if (config.obsPhi[pObsIndex] < 0.0)
        while (config.obsPhi[pObsIndex] < 0.0) config.obsPhi[pObsIndex] += (2.0 * PI);

      if (config.obsPhi[pObsIndex] > (2.0 * PI) )
        while (config.obsPhi[pObsIndex] > (2.0 * PI) ) config.obsPhi[pObsIndex] -= (2.0 * PI);

  }

  // reset the offset to zero
  azi_sun -= phiOffset;
  phiOffset = 0.0;

}
/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/

/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
/*--*/    Vec_t                                                     /*---*/
/*--*/    vMas(Vec_t r, Node_t node )                               /*---*/
/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
{

  SphVec_t rSphAu;
  Vec_t velocity;

  rSphAu = cartToSphPosAu(r);
  masGetNode(rSphAu, node);
  velocity = sphToCartVector(masNode.mhdV, r);

  return velocity;

}
/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/


/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
/*--*/    Vec_t                                                     /*---*/
/*--*/    rungeKuttaFlow( Vec_t r0, Scalar_t dt, Node_t node )      /*---*/
/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
{

  Vec_t k1, k2, k3, k4, v1, v2, v3;
  Vec_t r1;

  static Scalar_t sixth = 1.0 / 6.0;
  static Scalar_t two = 2.0;
  static Scalar_t half = 0.5;

  // find the forward displacement
  k1 = vectorScalarMult(vMas(r0, node), dt / config.rScale);
  v1 = vMas(vectorAddition(r0, vectorScalarMult(k1, half)), node);
  k2 = vectorScalarMult(v1, dt / config.rScale);
  v2 = vMas(vectorAddition(r0, vectorScalarMult(k2, half)), node);
  k3 = vectorScalarMult(v2, dt / config.rScale);
  v3 = vMas(vectorAddition(r0, k3), node);
  k4 = vectorScalarMult(v3, dt / config.rScale);

  r1.x = r0.x + sixth * (k1.x + two*k2.x + two*k3.x + k4.x);
  r1.y = r0.y + sixth * (k1.y + two*k2.y + two*k3.y + k4.y);
  r1.z = r0.z + sixth * (k1.z + two*k2.z + two*k3.z + k4.z);

  return r1;

}
/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/


/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
/*--*/    void                                                      /*---*/
/*--*/    masMoveNodes( Scalar_t dt )                               /*---*/
/*-----------------------------------------------------------------------*/
// The RK4 in here uses masGetNode - the previous step's s factor and file
// data are correct here since they have not been updated yet.
/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
{

  Index_t   face, row, col, shell, shell0, idx;
  Vec_t     r0, r1;
  Scalar_t  rmag;
  Node_t node;

  /* INNER_ACTIVE_SHELL on innermost proc left on inner boundary */
  /* on all other procs, the INNER_ACTIVE_SHELL moves forward    */

  if (mpi_rank == INNER_PROC) {
    shell0 = INNER_ACTIVE_SHELL + 1; }
  else {
    shell0 = INNER_ACTIVE_SHELL;}


  for (face = 0; face < NUM_FACES; face++ ){
    for (row = 0; row < FACE_ROWS; row++  ){
      for (col = 0; col < FACE_COLS; col++  ){
        for (shell = shell0; shell < LOCAL_NUM_SHELLS; shell++){

          idx = idx_frcs(face,row,col,shell);

          node = grid[idx];

          r0 = node.r;

          // store the current position as rOld
          grid[idx].rOld = node.r;

          // 4th order RungeKutta
          r1 = rungeKuttaFlow(r0, dt, node);

          rmag = sqrt( (r1.x*r1.x) + (r1.y*r1.y) + (r1.z*r1.z) );

          grid[idx].r = r1;
          grid[idx].rmag = rmag;

        }
      }
    }
  }

}
/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/


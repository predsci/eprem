/*-----------------------------------------------
 -- EMMREM: flow.c
 --
 -- Tools and routines related to the flow field.
 --
 -- ______________CHANGE HISTORY______________
 -- ___________________END CHANGE HISTORY_____________________
 ------------------------------------------------*/

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
#include "flow.h"
#include "simCore.h"
#include "geometry.h"
#include "readMAS.h"
#include "masInterp.h"


/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
/*--*/    void                                                      /*---*/
/*--*/    updateMhd()                                               /*---*/
/*--*                                                                *---*/
/*--* Calculate the Mhd quantities and put them in the data          *---*/
/*--* structure.                                                     *---*/
/*--*                                                                *---*/
/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
{

  Index_t face, row, col, shell, idealShockNode, idx;
  SphVec_t radpos;

  for (face  = 0; face < NUM_FACES; face++ )
  {
    for (row   = 0; row < FACE_ROWS; row++ )
    {
      for (col   = 0; col < FACE_COLS; col++ )
      {
        for (shell = INNER_ACTIVE_SHELL; shell < LOCAL_NUM_SHELLS; shell++ )
        {

          idx = idx_frcs(face,row,col,shell);

          // Save current time MHD values for use with derivatives below.
          grid[idx].mhdBmagOld       = grid[idx].mhdBmag;
          grid[idx].mhdDensityOld    = grid[idx].mhdDensity;
          grid[idx].mhdVsphOld.r     = grid[idx].mhdVr;
          grid[idx].mhdVsphOld.theta = grid[idx].mhdVtheta;
          grid[idx].mhdVsphOld.phi   = grid[idx].mhdVphi;
          grid[idx].dsOld            = grid[idx].ds;
          grid[idx].rOlder           = grid[idx].r;

          // Check for coupling to MAS
          if (mhdGridStatus == MHD_MAS) {

            // Get spherical position.
            radpos = cartToSphPosAu(grid[idx].r);

            // Interpolate MAS MHD values to current position at desired time (tGlobal+dt in main loop).
            // NOTE: This does NOT set values in grid, only in temp global masNode.
            masGetNode(radpos, grid[idx]);

          }

          // Check to see if node is in ideal shock domain
          idealShockNode = 0;
          if (config.idealShock > 0) {
            radpos = cartToSphPosAu(grid[idx].r);
            if (radpos.r <= config.idealShockSpeed * t_global) {
              if (config.idealShockWidth == 0.0)
                idealShockNode = 1;
              else {
                if ((fmin((2 * PI) - fabs(radpos.phi - config.idealShockPhi), fabs(radpos.phi - config.idealShockPhi)) <
                    (config.idealShockWidth / (double)2.0)) &&
                    (fabs(radpos.theta - config.idealShockTheta) < (config.idealShockWidth / 2.0)) )
                  idealShockNode = 1;
              }
            }
          }

          // Update Density.
          grid[idx].mhdDensity = mhdDensity(grid[idx].rmag, idealShockNode);

          // Update Div v.
          if (grid[idx].mhdDensity <= 0.0){
            grid[idx].mhdDivV=0.0;
          } else {
            grid[idx].mhdDivV=-1.0 * log(grid[idx].mhdDensity/grid[idx].mhdDensityOld)/config.tDel;
          }

          //update V
          mhdV(grid[idx].r,
               grid[idx].rmag,
               &grid[idx].mhdVr,
               &grid[idx].mhdVtheta,
               &grid[idx].mhdVphi,
               &grid[idx].mhdVmag,
               &grid[idx].mhdVvec,
               idealShockNode);

          // update B
          mhdB(grid[idx].r,
               grid[idx].rmag,
               &grid[idx].mhdBr,
               &grid[idx].mhdBtheta,
               &grid[idx].mhdBphi,
               &grid[idx].mhdBmag,
               &grid[idx].mhdBvec,
               grid[idx].mhdVr);

          // update del x B/B^2 (NOT implemented for helio yet!)
          if (config.useDrift > 0)
            mhdCurlBoverB2(grid[idx].r, &grid[idx].curlBoverB2,	grid[idx].mhdVr, idealShockNode);

          // calculate the time derivative terms
          if ((mpi_rank == 0) && (shell == INNER_ACTIVE_SHELL)) {

            grid[idx].mhdDuPar = 0.0;
            grid[idx].mhdDlnB = 0.0;
            grid[idx].mhdDlnN = 0.0;

          } else {

            grid[idx].mhdDuPar = (grid[idx].mhdVr     - grid[idx].mhdVsphOld.r)     * grid[idx].mhdBr
                                +(grid[idx].mhdVtheta - grid[idx].mhdVsphOld.theta) * grid[idx].mhdBtheta
                                +(grid[idx].mhdVphi   - grid[idx].mhdVsphOld.phi)   * grid[idx].mhdBphi;

            grid[idx].mhdDuPar /= grid[idx].mhdBmag;

            grid[idx].mhdDlnB  = log(grid[idx].mhdBmag/grid[idx].mhdBmagOld);

            grid[idx].mhdDlnN  = log(grid[idx].mhdDensity/grid[idx].mhdDensityOld);

          }
        }
      }
    }
  }
}
/*------ END  updateMhd ( ) ----------------------------------------*/
/*------------------------------------------------------------------*/


/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
/*--*/    Vec_t                                                     /*---*/
/*--*/    delrFlow( Vec_t position, Scalar_t rmag, Node_t node,     /*---*/
/*---*/             Scalar_t dt)                                    /*---*/
/*--*                                                                *---*/
/*--* Flow field pushing nodes around.                               *---*/
/*--*                                                                *---*/
/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
{
  Vec_t    delR;
  Scalar_t scale;

  if (mhdGridStatus == MHD_DEFAULT) {

    delR.x   = dt * node.mhdVvec.x / config.rScale ;
    delR.y   = dt * node.mhdVvec.y / config.rScale ;
    delR.z   = dt * node.mhdVvec.z / config.rScale ;

  } else if (mhdGridStatus == MHD_MAS) {

    delR.x   = dt * node.mhdVvec.x / config.rScale ;
    delR.y   = dt * node.mhdVvec.y / config.rScale ;
    delR.z   = dt * node.mhdVvec.z / config.rScale ;

  } else {

    scale = ( config.mhdUs / ( rmag * config.rScale ) );
    delR.x   = dt * scale * position.x;
    delR.y   = dt * scale * position.y;
    delR.z   = dt * scale * position.z;

  }

  return( delR );

} /*------ END  delrFlow ( ) ---------------------------------*/
/*-----------------------------------------------------------------------*/

/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
/*--*/    Scalar_t                                                  /*---*/
/*--*/    divV( Scalar_t rmag,                                      /*---*/
/*--*/          Index_t idealShockNode )                             /*--*/
/*--*                                                                *---*/
/*--* Returns Divergence of the Flow field.                          *---*/
/*--* We return a number in units of 1/TAU where TAU is              *---*/
/*--* is the light crossing time of 1 AU. Currently, position is     *---*/
/*--* not used. THIS FUNCTION IS NOT IN USE PRESENTLY.               *---*/
/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
{
  Scalar_t div;
  Scalar_t mhdUs;
  SphVec_t V;

  if (mhdGridStatus == MHD_DEFAULT) {
    div = 2.0 * config.mhdUs  / ( rmag * config.rScale );

    // ideal shock test
    if ( (config.idealShock > 0) && (idealShockNode > 0) )
      div *= idealShockFactor(rmag * config.rScale);

  }
  else if (mhdGridStatus == MHD_MAS) {
    V = masNode.mhdV;;
    mhdUs=sqrt(V.r*V.r+V.theta*V.theta+V.phi*V.phi);
    div = 2.0 * mhdUs/ ( rmag * config.rScale );
  }
  else {
    div = 2.0 * config.mhdUs / ( rmag * config.rScale );
  }

  return( div );

} /*------ END  divV ( ) ------------------------------------------------*/
/*-----------------------------------------------------------------------*/

/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
/*--*/    Scalar_t                                                  /*---*/
/*--*/    mhdDensity(Scalar_t rmag,                                 /*---*/
/*--*/               Index_t idealShockNode)                        /*---*/
/*--*                                                                *---*/
/*--* Calculate the local density.                                   *---*/
/*--*                                                                *---*/
/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
{
  Scalar_t dens, rr;

  rr = rmag * config.rScale;

  if (mhdGridStatus == MHD_DEFAULT) {

    dens = config.mhdNsAu / ( rr*rr );

    // ideal shock test
    if ( (config.idealShock > 0) && (idealShockNode > 0) )
      dens *= idealShockFactor(rr);

  }
  else if (mhdGridStatus == MHD_MAS) {

    dens = masNode.mhdD;

  }
  else {

    dens = config.mhdNsAu / ( rr*rr );

  }

  return( dens );

} /*------ END  mhdDensity ( ) ------------------------------------------*/
/*-----------------------------------------------------------------------*/

/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
/*--*/    void                                                      /*---*/
/*--*/    mhdB( Vec_t position, Scalar_t rmag,                      /*---*/
/*--*/	        Scalar_t *Br,   Scalar_t *Btheta,                   /*---*/
/*--*/          Scalar_t *Bphi, Scalar_t *Bmag,                     /*---*/
/*--*/ 	        Vec_t *Bvec,    Scalar_t Vr)                        /*---*/
/*--*                                                                *---*/
/*--* Calculate the local vec(B).                                    *---*/
/*--*                                                                *---*/
/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
{
  Scalar_t rr, theta, phi;
  SphVec_t B;

  rr = rmag * config.rScale;

  theta = acos(position.z / rmag);
  phi = atan2(position.y, position.x);

  if (mhdGridStatus == MHD_DEFAULT) {

    *Br = config.mhdBsAu / (rr*rr);

    *Btheta = 0.0;

    *Bphi = -1.0 * rr * (*Br) * (config.omegaSun / (Vr + VERYSMALL) ) * sin(theta);


    // Note on ideal shock: The normal component of the field is unchanged
    // across the shock. The tangential component compresses but does not
    // rotate. Since Btheta == 0.0, only Bphi changes. Across the shock,
    //
    // Vr1 Bphi1 = Vr2 Bphi2 <=> Bphi2 = Bphi1 (Vr1/Vr2)
    //
    // The given value of Vr will thus self-consistently scale Bphi to account
    // for shock affects (without the need for additional if..else logic), as
    // long as updateMhd calls this routine after calling mhdV.

    // ideal shock test
    /*
    if ( (config.idealShock > 0) && (idealShockNode > 0) ) {

      *Br     *= idealShockFactor(rr);
      *Btheta *= idealShockFactor(rr);
      *Bphi   *= idealShockFactor(rr);

    }
    */

  }
  else if (mhdGridStatus == MHD_MAS) {
    B = masNode.mhdB;
    *Br = B.r;
    *Btheta = B.theta;
    *Bphi = B.phi;
  }

  *Bmag = sqrt( (*Br) * (*Br) + (*Btheta) * (*Btheta) + (*Bphi) * (*Bphi) );

  (*Bvec).x = (*Br)     * sin(theta) * cos(phi)
  +       (*Btheta) * cos(theta) * cos(phi)
  -       (*Bphi)   *              sin(phi);

  (*Bvec).y = (*Br)     * sin(theta) * sin(phi)
  +       (*Btheta) * cos(theta) * sin(phi)
  +       (*Bphi)   *              cos(phi);

  (*Bvec).z = (*Br)     * cos(theta)
  -       (*Btheta) * sin(theta) ;

} /*------ END  mhdB ( ) ------------------------------------------------*/
/*-----------------------------------------------------------------------*/

/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
/*--*/    void                                                      /*---*/
/*--*/    mhdV( Vec_t position, Scalar_t  rmag,                     /*---*/
/*--*/	        Scalar_t *Vr,   Scalar_t *Vtheta,                   /*---*/
/*--*/	        Scalar_t *Vphi, Scalar_t *Vmag,                     /*---*/
/*--*/          Vec_t *Vvec, Index_t idealShockNode)                /*---*/
/*--*                                                                *---*/
/*--* Calculate the local vec(V) independently of position and rmag  *---*/
/*--* (currently).                                                   *---*/
/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
{
  SphVec_t V;
  Scalar_t theta, phi;

  theta = acos(position.z / rmag);
  phi = atan2(position.y, position.x);

  if (mhdGridStatus == MHD_DEFAULT) {

    *Vr = config.mhdUs;

    if ( (config.idealShock > 0) && (idealShockNode > 0) )
      *Vr *= idealShockFactor(rmag * config.rScale);

    *Vtheta = 0.0;

    *Vphi = 0.0;

  } else if (mhdGridStatus == MHD_MAS) {
    V = masNode.mhdV;
    *Vr = V.r;
    *Vtheta = V.theta;
    *Vphi = V.phi;
  }

  *Vmag = sqrt( (*Vr) * (*Vr) + (*Vtheta) * (*Vtheta) + (*Vphi) * (*Vphi) );

  (*Vvec).x = (*Vr)     * sin(theta) * cos(phi)
  +       (*Vtheta) * cos(theta) * cos(phi)
  -       (*Vphi)   *              sin(phi);

  (*Vvec).y = (*Vr)     * sin(theta) * sin(phi)
  +       (*Vtheta) * cos(theta) * sin(phi)
  +       (*Vphi)   *              cos(phi);

  (*Vvec).z = (*Vr)     * cos(theta)
  -       (*Vtheta) * sin(theta) ;

}
/*---------------- END  mhdV ( ) ----------------------------------------*/
/*-----------------------------------------------------------------------*/


/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
/*--*/    void                                                      /*---*/
/*--*/    mhdCurlBoverB2( Vec_t r,                                  /*---*/
/*--*/										SphVec_t *curl,														/*---*/
/*--*/										Scalar_t Vr,                              /*---*/
/*--*/                    Index_t idealShockNode)                   /*---*/
/*--*                                                                *---*/
/*--*			pull in the curl of B over B2                              *---*/
/*--*                                                                *---*/
/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
{

  if (mhdGridStatus == MHD_MAS)
  {

    *curl = masNode.curlBoverB2;

  }
  else
  {

    *curl = curlBoverB2(r, Vr, idealShockNode);

  }

}
/*---------------- END  mhdV ( ) ----------------------------------------*/
/*-----------------------------------------------------------------------*/


/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
/*--*/    SphVec_t                                                  /*---*/
/*--*/    parkerB(	Scalar_t r, Scalar_t theta,											/*---*/
/*--*/							Scalar_t Vr, Index_t idealShockNode )           /*---*/
/*--*/							                                                /*---*/
/*--*                                                                *---*/
/*--*			return the parker spiral solution for the magnetic field   *---*/
/*--*                                                                *---*/
/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
{

  SphVec_t B;


  B.r = config.mhdBsAu / (r * r);

  B.theta = 0.0;

  B.phi = -1.0 * r * B.r * (config.omegaSun / (Vr + VERYSMALL) ) * sin(theta);

  if ( (config.idealShock > 0) && (idealShockNode > 0) ) {

    B.r     *= idealShockFactor(r);
    B.theta *= idealShockFactor(r);
    B.phi   *= idealShockFactor(r);

  }

  return B;

}
/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/

/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
/*--*/    SphVec_t                                                  /*---*/
/*--*/    curlBoverB2(	Vec_t rCart,              									/*---*/
/*--*/									Scalar_t Vr,																/*---*/
/*--*/                  Index_t idealShockNode)											/*---*/
/*--*/							                                                /*---*/
/*--*                                                                *---*/
/*--*			return del x B/B^2 using the parkerB() solution            *---*/
/*--*                                                                *---*/
/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
{

  Scalar_t rCellWidthAU, thetaCellWidth, phiCellWidth;

  SphVec_t r, rMinus, rPlus, thetaMinus, thetaPlus, phiMinus, phiPlus;
  SphVec_t B_rm, B_rp, B_tm, B_tp, B_pm, B_pp;

  Scalar_t Bmag_rm, Bmag_rp, Bmag_tm, Bmag_tp, Bmag_pm, Bmag_pp;

  Scalar_t rAU, A, B;

  SphVec_t curl;

  // Set radial and angular offsets.
  //   NOTE -- SHOULD BE A CONFIGURATION OPTION OR AUTOMATED AND NOT STATIC
  rCellWidthAU = 0.001;
  thetaCellWidth = PI / 100.0;
  phiCellWidth = PI / 100.0;

  // Convert to spherical and in AU.
  r.r = config.rScale * sqrt( rCart.x * rCart.x + rCart.y * rCart.y + rCart.z * rCart.z );
  r.theta = acos( rCart.z * config.rScale / r.r );
  r.phi = atan2(rCart.y, rCart.x);

  // Assign the spatial values for each delta direction and take into account
  // the boundaries for phi and theta.
  rPlus.r = r.r + rCellWidthAU;
  rPlus.theta = r.theta;
  rPlus.phi = r.phi;

  rMinus.r = r.r - rCellWidthAU;
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
  phiPlus.phi = r.phi;
  phiPlus.phi += phiCellWidth;

  if ( phiPlus.phi > (2.0 * PI) )
    phiPlus.phi -= (2.0 * PI);

  phiMinus.r = r.r;
  phiMinus.theta = r.theta;
  phiMinus.phi = r.phi;
  phiMinus.phi -= phiCellWidth;

  if ( phiMinus.phi < 0.0 )
    phiMinus.phi += (2.0 * PI);


  // now that we have the spatial components, collect up the actual field values
  B_rm = parkerB(rMinus.r, rMinus.theta, Vr, idealShockNode);
  B_rp = parkerB(rPlus.r, rPlus.theta, Vr, idealShockNode);

  B_tm = parkerB(thetaMinus.r, thetaMinus.theta, Vr, idealShockNode);
  B_tp = parkerB(thetaPlus.r, thetaPlus.theta, Vr, idealShockNode);

  B_pm = parkerB(phiMinus.r, phiMinus.theta, Vr, idealShockNode);
  B_pp = parkerB(phiPlus.r, phiPlus.theta, Vr, idealShockNode);


  // calculate the magnitudes
  Bmag_rm = sqrt( B_rm.r * B_rm.r + B_rm.theta * B_rm.theta + B_rm.phi * B_rm.phi);
  Bmag_rp = sqrt( B_rp.r * B_rp.r + B_rp.theta * B_rp.theta + B_rp.phi * B_rp.phi);

  Bmag_tm = sqrt( B_tm.r * B_tm.r + B_tm.theta * B_tm.theta + B_tm.phi * B_tm.phi);
  Bmag_tp = sqrt( B_tp.r * B_tp.r + B_tp.theta * B_tp.theta + B_tp.phi * B_tp.phi);

  Bmag_pm = sqrt( B_pm.r * B_pm.r + B_pm.theta * B_pm.theta + B_pm.phi * B_pm.phi);
  Bmag_pp = sqrt( B_pp.r * B_pp.r + B_pp.theta * B_pp.theta + B_pp.phi * B_pp.phi);


  // finally!  calculating the curl itself
  rAU = r.r;

  // r
  A = (0.5 / thetaCellWidth) * ( B_tp.phi * sin(thetaPlus.theta) / (Bmag_tp * Bmag_tp) -
                                B_tm.phi * sin(thetaMinus.theta) / (Bmag_tm * Bmag_tm) );

  B = (0.5 / phiCellWidth) * ( B_pp.theta / (Bmag_pp * Bmag_pp) - B_pm.theta / (Bmag_pm * Bmag_pm) );

  curl.r = (1.0 / (rAU * sin(r.theta))) * (A - B);

  // theta
  A = (1.0 / sin(r.theta)) * (0.5 / phiCellWidth) * ( B_pp.r / (Bmag_pp * Bmag_pp) - B_pm.r / (Bmag_pm * Bmag_pm) );

  B = (0.5 / rCellWidthAU) * ( rPlus.r * B_rp.phi / (Bmag_rp * Bmag_rp) -
                              rMinus.r * B_rm.phi / (Bmag_rm * Bmag_rm) );

  curl.theta = (1.0 / rAU) * (A - B);

  // phi
  A = (0.5 / rCellWidthAU) * ( rPlus.r * B_rp.theta / (Bmag_rp * Bmag_rp) -
                              rMinus.r * B_rm.theta / (Bmag_rm * Bmag_rm) );

  B = (0.5 / thetaCellWidth) * ( B_tp.r / (Bmag_tp * Bmag_tp)  - B_tm.r / (Bmag_tm * Bmag_tm) );

  curl.phi = (1.0 / rAU) * (A - B);

  return curl;

}
/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/

/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
/*--*/    Scalar_t                                                  /*---*/
/*--*/    idealShockFactor( Scalar_t r )                            /*---*/
/*--*                                                                *---*/
/*--*     Multiplicative factor for ideal shock fields (r in AU)     *---*/
/*--*                                                                *---*/
/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
{

  Scalar_t factor, exponent, ratio, falloff;

  exponent = exp(2.0 *
                 (config.idealShockSharpness / config.idealShockScaleLength) *
                 (r - ((config.idealShockSpeed / C) * (t_global - config.idealShockInitTime / DAY) + config.rScale)));

  if (exponent >= 1.0e8)
    ratio = 1.0;
  else if (exponent <= 1.0e-8)
    ratio = -1.0;
  else
    ratio = (exponent - 1.0) / (exponent + 1.0);

  factor = (1.0 + ((1.0 - ratio) * 0.5 * (config.idealShockJump - 1.0)));
  falloff = exp(config.idealShockFalloff *
    (r - ((config.idealShockSpeed / C) * (t_global - config.idealShockInitTime / DAY) + config.rScale))
  );

  return 1 + (factor - 1) * falloff;

}
/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/

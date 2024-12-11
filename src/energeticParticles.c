/*-----------------------------------------------
 -- EMMREM: energeticParticles.c
 --
 -- updateEnergeticParticles():
 --     Functionality to update energetic particle distributions.
 --
 -- ______________CHANGE HISTORY______________
 --
 -- ______________END CHANGE HISTORY______________
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
#include <stdio.h>
#include <float.h>
#include <stdlib.h>
#include <string.h>

#include "baseTypes.h"
#include "global.h"
#include "configuration.h"
#include "energeticParticles.h"
#include "energeticParticlesBoundary.h"
#include "unifiedOutput.h"
#include "simCore.h"
#include "geometry.h"
#include "error.h"
#include "flow.h"
#include "readMAS.h"
#include "masInterp.h"
#include "observerOutput.h"
#include "timers.h"

Scalar_t *deltaShell;
Scalar_t *shockDist;

Index_t streamlistSize;
Index_t *shockFlag;
Index_t *shockCalcFlag;

Scalar_t dsMin;

Index_t maxsubcycles_energychange = 0;
Index_t maxsubcycles_focusing = 0;
Index_t maxsubcycles_energychangeGlobal = 0;
Index_t maxsubcycles_focusingGlobal = 0;
Scalar_t min_tau        = DBL_MAX;
Scalar_t min_tau_global = DBL_MAX;

Scalar_t leaving_left = 0;
Scalar_t leaving_right = 0;
Scalar_t leaving_leftGlobal = 0;
Scalar_t leaving_rightGlobal = 0;

/*---------------------------------------------------------------*/
/*---------------------------------------------------------------*/
/*--*/    Scalar_t                                           /*--*/
/*--*/    meanFreePath(Index_t species,                      /*--*/
/*--*/                 Index_t energy,                       /*--*/
/*--*/                 Scalar_t range)                       /*--*/
/*--*/                                                       /*--*/
/*--  range needs to be expressed in AU, not code units        --*/
/* -- Mean Free Path in AU ------------------------------------- */
{/*--------------------------------------------------------------*/

  Scalar_t mfp;

  mfp = rigidity[idx_se(species, energy)]
        * pow(range, config.mfpRadialPower) * config.lamo;

  return mfp;

} /*----------------------------------------------------------------*/
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/*--*/             void                                         /*--*/
/*--*/     updateEnergeticParticles( void )                     /*--*/
/*--                                                              --*/
/*-- Evaluate current particle distributions for all energy       --*/
/*-- ranges, for all shells and every shell node in each.         --*/
/*-- Called by the main simulation loop every time step, just     --*/
/*-- before the shell is propagated in space by the flow field.   --*/
/*-- NB-The INNER_SHELL does not need to be included in the       --*/
/*-- update as its data will be overwritten when it is used       --*/
/*-- as a template to spawn a new shell.                          --*/
/*------------------------------------------------------------------*/
{/*-----------------------------------------------------------------*/

  Index_t face, row, col, shell, innerComputeShell, step, idx;
  Time_t  t_global_saved;
  Scalar_t dt,tau;
  Index_t computeIndex, numIters, iterIndex, species, energy;

  // Each MPI rank computes numIters number of streams seqentially.
  numIters = NUM_STREAMS / N_PROCS;

  double timer_tmp = 0;

  /* Allocate temporary arrays used in several functions here */
  deltaShell    = (Scalar_t*) malloc(sizeof(Scalar_t)*FRC*NUM_MUSTEPS);
  shockDist     = (Scalar_t*) malloc(sizeof(Scalar_t)*NUM_FACES
                                                  *FACE_ROWS*FACE_COLS
                                         *LOCAL_NUM_SHELLS*NUM_SPECIES
                                               *NUM_ESTEPS*NUM_MUSTEPS);
  shockFlag     = (Index_t*)  calloc(NUM_FACES*FACE_ROWS*FACE_COLS
                                     *LOCAL_NUM_SHELLS*NUM_SPECIES,
                                     sizeof(Index_t));
  shockCalcFlag = (Index_t*)  calloc(NUM_FACES*FACE_ROWS*FACE_COLS
                                     *LOCAL_NUM_SHELLS*NUM_SPECIES,
                                     sizeof(Index_t));

  // Save global time (the time step update happens after this routine).
  t_global_saved = t_global;

  // Calculate the sub-timestep to use in the mhd and node movement.
  dt = config.tDel / (1.0 * config.numEpSteps);

  // Loop over the number of EP steps.
  for (step = 0; step < config.numEpSteps; step++ )
  {

    // Requires entire stream on one process.  Sequential on the rank.
    for (iterIndex = 0; iterIndex <= numIters; iterIndex++)
    {

      // Gather up current stream from shells across all ranks.
      update_stream_from_shells( iterIndex );

      // Set stream values that require +/- nodes along the stream.
      // NOTE!  This sets imporant values used in focusing!
      updateStreamValues( iterIndex );

      if (config.checkSeedPopulation > 0){
        seedPopulationCheck( iterIndex );
      }

      if (config.useParallelDiffusion > 0)
      {
        timer_tmp = MPI_Wtime();

        // Thin out streams based on dsh parameters
        // (ignore streams too close together).
        GetStreamList( iterIndex );

        // Compute the stream "diffusion" (advection along the stream).
        DiffuseStreamData( iterIndex, dt );

        timer_diffusestream = timer_diffusestream
                              + (MPI_Wtime() - timer_tmp);
      }

      // Scatter current stream to shells across all ranks.
      update_shells_from_stream( iterIndex );

    }

    // Make sure not to compute the inner shell on proc 0
    // since it has no time history
    if (mpi_rank == 0) {
      innerComputeShell = INNER_ACTIVE_SHELL + 1;
    } else {
      innerComputeShell = INNER_ACTIVE_SHELL;
    }

    for (shell = innerComputeShell; shell < LOCAL_NUM_SHELLS; shell++ )
    {

      // computed on local process
      for (computeIndex = 0; computeIndex < NUM_STREAMS; computeIndex++)
      {

        face = computeLines[computeIndex][0];
        row  = computeLines[computeIndex][1];
        col  = computeLines[computeIndex][2];
        idx = idx_frcs(face,row,col,shell);
//
//     ****** Find minimum mean free path time scale.
//
        for (species = 0; species < NUM_SPECIES; species++) {
          for (energy = 0; energy < NUM_ESTEPS; energy++) {
            tau = meanFreePath(species, energy,
                  grid[idx_frcs(face,row,col,shell)].rmag*config.rScale)
                  /vgrid[energy];
            if (tau < min_tau){
              min_tau = tau;
            }
          }
        }
//
//      ****** ADIABATIC FOCUS ******
//
        if ( config.useAdiabaticFocus > 0 ){

          timer_tmp = MPI_Wtime();

          AdiabaticFocusing(face, row, col, shell, dt);

          timer_adiabaticfocus = timer_adiabaticfocus
                                 + (MPI_Wtime() - timer_tmp);

        }
//
//      ****** ADIABATIC CHANGE ******
//
        if ( config.useAdiabaticChange > 0 ) {

          timer_tmp = MPI_Wtime();

          // If not running MAS-coupled run, always accelerate
          // everywhere in the domain.
          if (config.masCouple == 0) {
            AdiabaticChange(face, row, col, shell, dt);
          }
          //MAS Coupling Run.
          else {
            if (config.masHelCouple == 0){
              AdiabaticChange(face, row, col, shell, dt);
            }
            // If this is a heliospheric run, only accelerate
            // in corona when it is still evolving.
            else {
              if (t_global <= masTime[config.masNumFiles - 1]){
                AdiabaticChange(face, row, col, shell, dt);
              }
              else if ( (grid[idx].rmag * config.rScale)
                      > (config.masRadialMax * RSAU)     ){
                AdiabaticChange(face, row, col, shell, dt);
              }
            }
          }

          timer_adiabaticchange = timer_adiabaticchange
                                  + (MPI_Wtime() - timer_tmp);

        } // adiabaticChange

      } // Stream index
//
//    ****** DIFFUSE SHELL DATA ******
//
      if ( config.useShellDiffusion > 0){

        timer_tmp = MPI_Wtime();

        DiffuseShellData( shell, dt );

        timer_diffuseshell = timer_diffuseshell
                           + (MPI_Wtime() - timer_tmp);
      }
//
//    ****** DRIFT SHELL DATA ******
//
      if (config.useDrift > 0){

        timer_tmp = MPI_Wtime();

        DriftShellData( shell, dt );

        timer_driftshell = timer_driftshell
                           + (MPI_Wtime() - timer_tmp);

      }

    }
//
//   ****** UPDATE EpSubcycle TIME ******
//
    t_global += dt;
  }

  // Free up temporary arrays.
  free(deltaShell);
  free(shockDist);
  free(shockFlag);
  free(shockCalcFlag);

  // Reset time since t_global is updated after this routine.
  t_global = t_global_saved;

}/*-------- END updateEnergeticParticles() -------------------------*/
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/*--*/  Vec_t                                                   /*--*/
/*--*/  driftVelocity( Index_t species,                         /*--*/
/*--*/                 Index_t energy,                          /*--*/
/*--*/                 SphVec_t curlBoverB2,                    /*--*/
/*--*/                 Vec_t r )                                /*--*/
/*--*/                                                          /*--*/
/*--*/                                                          /*--*/
/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
{

  Scalar_t p, factor, om1;

  SphVec_t vd;

  Vec_t vd_cart;

  p = exp(lnpmin + energy * dlnp);

  om1 = config.charge[species] * OM / config.mass[species];

  factor = (1.0 / om1) * (p * vgrid[energy] / 3.0);

/*   OM is e B0 / m c and is normalized by AU/C  */
/*   expects curlBoverB2 in units [ (1 / MHD_B_NORM ) x  ( 1 / AU ) ] */
/*   in other words curlBoverB2 should use all code unit fields   */
/*   and gradient in 1/AU. So in curlBoverB2, need to normalize distances by config.Rscale */
/*   but do not normalize magnetic field by MHD_B_NORM     */

  vd.r = factor * curlBoverB2.r;
  vd.theta = factor * curlBoverB2.theta;
  vd.phi = factor * curlBoverB2.phi;

  vd_cart = sphToCartVector(vd, r);

  return vd_cart;

}
/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/*--*/    void                                                  /*--*/
/*--*/    ShellData( void )                                     /*--*/
/*--                                                              --*/
/*-- The perpendicular lengths to NEWS neighbors                  --*/
/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
{

  Index_t shell, face, row, col;
  Vec_t r, r1, en;
  Vec_t eb;

  Node_t     node, node1;
  Neighbor_t n, e, w, s;

  Scalar_t dl;
  Scalar_t dot;

  for (shell = INNER_ACTIVE_SHELL; shell < LOCAL_NUM_SHELLS; shell++)
  {

    dlPerMin[shell] = 1.0e20;

    for (face = 0; face < NUM_FACES; face++) {
      for (row  = 0; row < FACE_ROWS; row++) {
        for (col  = 0;  col < FACE_COLS; col++) {

          node = grid[idx_frcs(face,row,col,shell)];
          n    = node.n;
          e    = node.e;
          w    = node.w;
          s    = node.s;

          r = node.r;

          r.x *= config.rScale;
          r.y *= config.rScale;
          r.z *= config.rScale;

          eb.x = node.mhdBvec.x / node.mhdBmag;
          eb.y = node.mhdBvec.y / node.mhdBmag;
          eb.z = node.mhdBvec.z / node.mhdBmag;

          /* North Neighbor */
          node1 = grid[idx_frcs(n.face,n.row,n.col,n.shell)];
          r1 = node1.r;

          r1.x *= config.rScale;
          r1.y *= config.rScale;
          r1.z *= config.rScale;

          en.x = r1.x - r.x;
          en.y = r1.y - r.y;
          en.z = r1.z - r.z;

          dl = sqrt(en.x*en.x + en.y*en.y + en.z*en.z);
          grid[idx_frcs(face,row,col,shell)].n.dl = dl;

          en.x /= dl;
          en.y /= dl;
          en.z /= dl;

          dot = en.x*eb.x + en.y*eb.y + en.z*eb.z;

          grid[idx_frcs(face,row,col,shell)].n.dlPer =
          dl * sqrt( 1.0 - dot * dot );


          if (grid[idx_frcs(face,row,col,shell)].n.dlPer < dlPerMin[shell]) {
            dlPerMin[shell] = grid[idx_frcs(face,row,col,shell)].n.dlPer;
          }

          /* South neighbor */
          node1 = grid[idx_frcs(s.face,s.row,s.col,s.shell)];
          r1 = node1.r;

          r1.x *= config.rScale;
          r1.y *= config.rScale;
          r1.z *= config.rScale;

          en.x = r1.x - r.x;
          en.y = r1.y - r.y;
          en.z = r1.z - r.z;

          dl = sqrt(en.x*en.x + en.y*en.y + en.z*en.z);
          grid[idx_frcs(face,row,col,shell)].s.dl = dl;

          en.x /= dl;
          en.y /= dl;
          en.z /= dl;

          dot = en.x*eb.x + en.y*eb.y + en.z*eb.z;

          grid[idx_frcs(face,row,col,shell)].s.dlPer =
          dl * sqrt( 1.0 - dot * dot );

          if (grid[idx_frcs(face,row,col,shell)].s.dlPer  < dlPerMin[shell]) {
            dlPerMin[shell] = grid[idx_frcs(face,row,col,shell)].s.dlPer;
          }

          /* E Neighbor */
          node1 = grid[idx_frcs(e.face,e.row,e.col,e.shell)];
          r1 = node1.r;

          r1.x *= config.rScale;
          r1.y *= config.rScale;
          r1.z *= config.rScale;

          en.x = r1.x - r.x;
          en.y = r1.y - r.y;
          en.z = r1.z - r.z;

          dl = sqrt(en.x*en.x + en.y*en.y + en.z*en.z);
          grid[idx_frcs(face,row,col,shell)].e.dl = dl;

          en.x /= dl;
          en.y /= dl;
          en.z /= dl;

          dot = en.x*eb.x + en.y*eb.y + en.z*eb.z;

          grid[idx_frcs(face,row,col,shell)].e.dlPer =
          dl * sqrt( 1.0 - dot * dot );

          if (grid[idx_frcs(face,row,col,shell)].e.dlPer  < dlPerMin[shell]) {
            dlPerMin[shell] = grid[idx_frcs(face,row,col,shell)].e.dlPer;
          }

          /* W Neighbor */
          node1 = grid[idx_frcs(w.face,w.row,w.col,w.shell)];
          r1 = node1.r;

          r1.x *= config.rScale;
          r1.y *= config.rScale;
          r1.z *= config.rScale;

          en.x = r1.x - r.x;
          en.y = r1.y - r.y;
          en.z = r1.z - r.z;

          dl = sqrt(en.x*en.x + en.y*en.y + en.z*en.z);
          grid[idx_frcs(face,row,col,shell)].w.dl = dl;

          en.x /= dl;
          en.y /= dl;
          en.z /= dl;

          dot = en.x*eb.x + en.y*eb.y + en.z*eb.z;

          grid[idx_frcs(face,row,col,shell)].w.dlPer = dl * sqrt( 1.0 - dot * dot );

          if (grid[idx_frcs(face,row,col,shell)].w.dlPer  < dlPerMin[shell]) {
            dlPerMin[shell] = grid[idx_frcs(face,row,col,shell)].w.dlPer;
          }

        }
      }
    }
  }

}/*-------- END ShellData() ----------------------------------------*/
/*------------------------------------------------------------------*/

/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/*--*/     void                                                 /*--*/
/*--*/     DiffuseShellData(  Index_t shell,                    /*--*/
/*--*/                        Scalar_t dt )                     /*--*/
/*--                                                              --*/
/*-- Diffuse Data within a Shell                                  --*/
/*-- perp diffusion                                               --*/
/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
{/*-----------------------------------------------------------------*/

  Index_t face, row, col, species, energy, mu;

  Scalar_t* mfp;

  Scalar_t dt_kper, delN, delE, delW, delS;

  Node_t node;

  mfp = (Scalar_t*)malloc(NUM_FACES*FACE_ROWS*FACE_COLS*sizeof(Scalar_t));

  for (species = 0; species < NUM_SPECIES; species++)
  {
    for (energy = 0; energy < NUM_ESTEPS; energy++)
    {

      // set mfp and zero out deltaShell
      for (face = 0; face < NUM_FACES; face++)
      {
        for (row = 0; row < FACE_ROWS; row++)
        {
          for (col = 0; col < FACE_COLS; col++)
          {

            mfp[idx_frc(face,row,col)] = 0.0;

            for (mu = 0; mu < NUM_MUSTEPS; mu++)
              deltaShell[idx_frcm(face,row,col,mu)] = 0.0;

            mfp[idx_frc(face,row,col)] =
               meanFreePath(species, energy,
               grid[idx_frcs(face,row,col,shell)].rmag * config.rScale);

          }

        }

      }

      // calculate over all faces, including observer faces
      // the neighbor connections internal to the observer faces remain fixed and
      //  only the edges of each observer face connect to "real" nodes
      for (face = 0; face < NUM_FACES; face++ )
      {
        for (row  = 0; row  < FACE_ROWS; row++ )
        {
          for (col = 0; col < FACE_COLS; col++ )
          {

            node = grid[idx_frcs(face,row,col,shell)];

            dt_kper = dt * 0.33333333 * mfp[idx_frc(face,row,col)]
              * vgrid[energy] * config.kperxkpar;

            delN = dt_kper / (node.n.dlPer * node.n.dlPer + VERYSMALL);
            if (delN > THRESH)
              delN = THRESH;

            delE = dt_kper / (node.e.dlPer * node.e.dlPer + VERYSMALL);
            if (delE > THRESH)
              delE = THRESH;

            delW = dt_kper / (node.w.dlPer * node.w.dlPer + VERYSMALL);
            if (delW > THRESH)
              delW = THRESH;

            delS = dt_kper / (node.s.dlPer * node.s.dlPer + VERYSMALL);
            if (delS > THRESH)
              delS = THRESH;

            for (mu = 0; mu < NUM_MUSTEPS; mu++)
            {

              deltaShell[idx_frcm(node.n.face,node.n.row,node.n.col,mu)]
              += delN * eParts[idx_frcsspem(face,row,col,shell,species,energy,mu)];

              deltaShell[idx_frcm(node.e.face,node.e.row,node.e.col,mu)]
              += delE * eParts[idx_frcsspem(face,row,col,shell,species,energy,mu)];

              deltaShell[idx_frcm(node.w.face,node.w.row,node.w.col,mu)]
              += delW * eParts[idx_frcsspem(face,row,col,shell,species,energy,mu)];

              deltaShell[idx_frcm(node.s.face,node.s.row,node.s.col,mu)]
              += delS * eParts[idx_frcsspem(face,row,col,shell,species,energy,mu)];

              deltaShell[idx_frcm(face,row,col,mu)]
              -= (delN+delE+delW+delS) * eParts[idx_frcsspem(face,row,col,shell,species,energy,mu)];

            }

          }

        }

      }


      for (face = 0; face < NUM_FACES; face++){
        for (row= 0; row < FACE_ROWS; row++){
          for (col = 0; col < FACE_COLS; col++){
            for (mu = 0; mu < NUM_MUSTEPS; mu++){

              eParts[idx_frcsspem(face,row,col,shell,species,energy,mu)]
                               += deltaShell[idx_frcm(face,row,col,mu)];

              // check for NaN and Inf
              // check for NaN and Inf
              checkNaN(mpi_rank, face, row, col, shell,
                       eParts[idx_frcsspem(face,row,col,shell,species,energy,mu)],
                       "DiffuseShellData");
              checkInf(mpi_rank, face, row, col, shell,
                       eParts[idx_frcsspem(face,row,col,shell,species,energy,mu)],
                       "DiffuseShellData");

            }
          }
        }
      }

    }

  }

  free(mfp);

}
/*--------- END DiffuseShellData()    ------------------------------*/
/*------------------------------------------------------------------*/

/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/*--*/     void                                                 /*--*/
/*--*/     DriftShellData(  Index_t shell,                      /*--*/
/*--*/              Scalar_t dt )                               /*--*/
/*--                                                              --*/
/*-- Move Data within a Shell                                     --*/
/*-- due to drift motion                                          --*/
/*------------------------------------------------------------------*/
/*-----------------------------------------------------------------*/
{

  Index_t species, energy, face, row, col, mu;

  Node_t node, node1;
  Neighbor_t n, e, w, s;

  Vec_t rv, rv1, en, vd;

  Scalar_t vdp, del;


  for (species = 0; species < NUM_SPECIES; species++)
  {
    for (energy = 0; energy < NUM_ESTEPS; energy++)
    {

      for (face = 0; face < NUM_FACES; face++ )
      {
        for (row = 0; row < FACE_ROWS; row++ )
        {
          for (col = 0; col < FACE_COLS; col++ )
          {
            for (mu = 0; mu < NUM_MUSTEPS; mu++ )
            {

              deltaShell[idx_frcm(face,row,col,mu)] = 0.0;

            }
          }
        }
      }

      for (face = 0; face < NUM_FACES; face++ )
      {
        for (row = 0; row < FACE_ROWS; row++ )
        {
          for (col = 0; col < FACE_COLS; col++ )
          {

            node = grid[idx_frcs(face,row,col,shell)];

            n    = node.n;
            e    = node.e;
            w    = node.w;
            s    = node.s;

            rv = node.r;

            vd = driftVelocity(species, energy, node.curlBoverB2, rv);

            /* N move */

            node1 = grid[idx_frcs(n.face,n.row,n.col,n.shell)];
            rv1 = node1.r;

            en.x = rv1.x - rv.x;
            en.y = rv1.y - rv.y;
            en.z = rv1.z - rv.z;

            vdp = (vd.x * en.x + vd.y * en.y + vd.z * en.z ) * config.rScale / node.n.dl;

            if (vdp > 0.0)
            {

              del = vdp * dt /  node.n.dl  ;

              if (del > THRESH)
                del = THRESH;

              //vd.x -= vdp * en.x / node.n.dl;
              //vd.y -= vdp * en.y / node.n.dl;
              //vd.z -= vdp * en.z / node.n.dl;

              for (mu = 0; mu < NUM_MUSTEPS; mu++)
              {

                deltaShell[idx_frcm(node.n.face,node.n.row,node.n.col,mu)]
                  += del * eParts[idx_frcsspem(face,row,col,shell,species,energy,mu)];

                deltaShell[idx_frcm(face,row,col,mu)]
                  -= del * eParts[idx_frcsspem(face,row,col,shell,species,energy,mu)];

              }

            }

            /* E move */

            node1 = grid[idx_frcs(e.face,e.row,e.col,e.shell)];
            rv1 = node1.r;

            en.x = rv1.x - rv.x;
            en.y = rv1.y - rv.y;
            en.z = rv1.z - rv.z;

            vdp = (vd.x * en.x + vd.y * en.y + vd.z * en.z ) * config.rScale / node.e.dl;

            if (vdp > 0.0)
            {

              del = vdp * dt /  node.e.dl  ;

              if (del > THRESH)
                del = THRESH;

              //vd.x -= vdp * en.x / node.e.dl;
              //vd.y -= vdp * en.y / node.e.dl;
              //vd.z -= vdp * en.z / node.e.dl;

              for (mu = 0; mu<NUM_MUSTEPS; mu++)
              {

                deltaShell[idx_frcm(node.e.face,node.e.row,node.e.col,mu)]
                  += del * eParts[idx_frcsspem(face,row,col,shell,species,energy,mu)];

                deltaShell[idx_frcm(face,row,col,mu)]
                  -= del * eParts[idx_frcsspem(face,row,col,shell,species,energy,mu)];

              }

            }

            /* W diffuse */

            node1 = grid[idx_frcs(w.face,w.row,w.col,w.shell)];
            rv1 = node1.r;

            en.x = rv1.x - rv.x;
            en.y = rv1.y - rv.y;
            en.z = rv1.z - rv.z;

            vdp = (vd.x * en.x + vd.y * en.y + vd.z * en.z ) * config.rScale / node.w.dl;

            if (vdp > 0.0)
            {

              del = vdp * dt /  node.w.dl ;

              if (del > THRESH)
                del = THRESH;

              //vd.x -= vdp * en.x / node.w.dl;
              //vd.y -= vdp * en.y / node.w.dl;
              //vd.z -= vdp * en.z / node.w.dl;

              for (mu = 0; mu<NUM_MUSTEPS; mu++)
              {

                deltaShell[idx_frcm(node.w.face,node.w.row,node.w.col,mu)]
                  += del * eParts[idx_frcsspem(face,row,col,shell,species,energy,mu)];

                deltaShell[idx_frcm(face,row,col,mu)]
                  -= del * eParts[idx_frcsspem(face,row,col,shell,species,energy,mu)];

              }

            }

            /* S diffuse */

            node1 = grid[idx_frcs(s.face,s.row,s.col,s.shell)];
            rv1 = node1.r;

            en.x = rv1.x - rv.x;
            en.y = rv1.y - rv.y;
            en.z = rv1.z - rv.z;

            vdp = (vd.x * en.x + vd.y * en.y + vd.z * en.z ) * config.rScale  / node.s.dl;

            if (vdp > 0.0)
            {

              del = vdp * dt /  node.s.dl ;

              if (del > THRESH)
                del = THRESH;

              //vd.x -= vdp*en.x/node.s.dl;
              //vd.y -= vdp*en.y/node.s.dl;
              //vd.z -= vdp*en.z/node.s.dl;

              for (mu = 0; mu<NUM_MUSTEPS; mu++)
              {

                deltaShell[idx_frcm(node.s.face,node.s.row,node.s.col,mu)]
                  += del * eParts[idx_frcsspem(face,row,col,shell,species,energy,mu)];

                deltaShell[idx_frcm(face,row,col,mu)]
                  -= del * eParts[idx_frcsspem(face,row,col,shell,species,energy,mu)];

              }

            }


          }

        }

      }


      for (face = 0; face < NUM_FACES; face++ )
      {
        for (row = 0; row < FACE_ROWS; row++ )
        {
          for (col = 0; col < FACE_COLS; col++ )
          {
            for (mu = 0; mu < NUM_MUSTEPS; mu++)
            {

              eParts[idx_frcsspem(face,row,col,shell,species,energy,mu)]
                += deltaShell[idx_frcm(face,row,col,mu)];

              // check for NaN and Inf
              checkNaN(mpi_rank, face, row, col, shell,
                       eParts[idx_frcsspem(face,row,col,shell,species,energy,mu)],
                       "DriftShellData");
              checkInf(mpi_rank, face, row, col, shell,
                       eParts[idx_frcsspem(face,row,col,shell,species,energy,mu)],
                       "DriftShellData");

            }
          }
        }
      }

    }
  }

}
/*----------- END DriftShellData()    ------------------------------*/
/*------------------------------------------------------------------*/

/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/*--*/             void                                         /*--*/
/*--*/     AdiabaticChange(  Index_t face,                      /*--*/
/*--*/                       Index_t row,                       /*--*/
/*--*/                       Index_t col,                       /*--*/
/*--*/                       Index_t shell,                     /*--*/
/*--*/                       Scalar_t dt )                      /*--*/
/*--                                                              --*/
/*-- Evaluate the adiabatic change using upwinding                --*/
/*------------------------------------------------------------------*/
{/*-----------------------------------------------------------------*/

// Want to solve f^n+1 = f^n - dt * v * df/dx
// Here, all d/dt MHD derivative terms have not yet been divided by dt.
// The "dt" passed into this function is the EP subcyled time step.
// Therefore, we calcualate the large DT (the full EPREM time-step) to
// compute the MHD derivative constants.

  Node_t node;

  Index_t N_subcycles                      = 1;
  Scalar_t safety_factor                   = 0.9;
  Scalar_t vel_abs_max                     = 0.0;
  Scalar_t half                            = 0.5;
  Scalar_t one                             = 1.0;
  Index_t s, species, energy, mu, idx;

  Scalar_t *restrict vel, *restrict f, *restrict f1;
  Scalar_t *restrict v_avg;

  Scalar_t dt_full, dt_stable, dt_subcycle;
  Scalar_t DlnBDt, DlnNDt, DuParDt;
  Scalar_t a, b, muval, cool_damp_fac;

  // Allocate space for the dist function, effective advection velocity, and fluxes.

  f    = (Scalar_t *) malloc(sizeof(Scalar_t) * NUM_SPECIES * NUM_ESTEPS * NUM_MUSTEPS);
  f1   = (Scalar_t *) malloc(sizeof(Scalar_t) * NUM_SPECIES * NUM_ESTEPS * NUM_MUSTEPS);
  vel  = (Scalar_t *) malloc(sizeof(Scalar_t) * NUM_SPECIES * NUM_ESTEPS * NUM_MUSTEPS);
  v_avg  = (Scalar_t *) malloc(sizeof(Scalar_t) * NUM_SPECIES * (NUM_ESTEPS+1) * NUM_MUSTEPS);

  // Get the current node.  This contains the MHD differences
  // after the node has been moved (e.g. Delta-MHD = MHD^n+1-MHD^n)

  node = grid[idx_frcs(face,row,col,shell)];

  // Get the full timestep value and use it to compute MHD derivative terms.
  dt_full = dt*config.numEpSteps;
  DlnBDt  = node.mhdDlnB/dt_full;
  DlnNDt  = node.mhdDlnN/dt_full;

  // Cooling damping:
  if (DlnNDt<0 && config.dampCooling>0){
    cool_damp_fac = config.dampCooling * 3.0/5.0;
  }else{
    cool_damp_fac = 1.0;
  }

  DuParDt = node.mhdDuPar/dt_full;

  // For accuracy, we advect ln(distribution), so need to convert here.
  // In same loop, we compute the effective advection velocity (ln(p)/time)
  // across the energy grid.
  // We keep track of the maximum |velocity| to compute the stable timestep.

  for (species = 0; species < NUM_SPECIES; species++) {
    for (energy = 0; energy < NUM_ESTEPS; energy++) {
      for (mu = 0; mu < NUM_MUSTEPS; mu++) {

        idx = idx_spem(species,energy,mu);

        f[idx] = log(eParts[idx_frcsspem(face,row,col,shell,species,energy,mu)]);

        muval = mugrid[mu];

        a = -muval * DuParDt;
        b = muval*muval * (DlnNDt - DlnBDt)
            + half*(one-muval*muval) * DlnBDt;

        vel[idx] = cool_damp_fac * (a/vgrid[energy] + b);

        if (fabs(vel[idx]) > vel_abs_max) vel_abs_max = fabs(vel[idx]);

      }
    }
  }

  // Generate extended v_avg array.
  // Eventually, and alternative vgrid array should be made to
  // match the "half mesh" directly so no averaging is needed.
  // Here, the CFL could be relaxed a bit if we used v_avg to set it
  // instead of vel above.

  for (species = 0; species < NUM_SPECIES; species++) {
    for (energy = 0; energy < NUM_ESTEPS-1; energy++) {
      for (mu = 0; mu < NUM_MUSTEPS; mu++) {

        v_avg[idx_spep1m(species,energy+1,mu)] = half * (
                          vel[idx_spem(species,energy+1,mu)]
                        + vel[idx_spem(species,energy  ,mu)]);
      }
    }
  }

  // Make extended v_avg boundary values (just copy for now).

  for (species = 0; species < NUM_SPECIES; species++) {
    for (mu = 0; mu < NUM_MUSTEPS; mu++) {
      energy=0;
      v_avg[idx_spep1m(species,energy,mu)] = v_avg[idx_spep1m(species,energy+1,mu)];
      energy=NUM_ESTEPS;
      v_avg[idx_spep1m(species,energy,mu)] = v_avg[idx_spep1m(species,energy-1,mu)];
    }
  }

  // Find the stable time-step based on the velocity and dlnp.
  // Here, we are relying on the fact that dlnp is NOT dependent on species
  // anymore, so we use only species index 0.

  dt_stable = safety_factor * dlnp/vel_abs_max;

  // Now that we have the stable timestep, we need to be careful
  // about how many subcycles to use.
  // This routine is being called within the
  // EpSubCycle loop, so the dt we need to go is dt_full/config.numEpSteps,
  // or, simply the "dt" sent into the routine.

  N_subcycles = (Index_t) ceil(dt/dt_stable);

  // Keep track of the maximum number of subcycles
  if (N_subcycles > maxsubcycles_energychange){
    maxsubcycles_energychange = N_subcycles;
  }

  // Set the dt_subcycle so that it exactly reaches dt.

  dt_subcycle = dt/N_subcycles;

  // Now compute the sub-cycled upwind advance.
  if (AdiabaticChangeAlg == 2 || AdiabaticChangeAlg == 3){

    Scalar_t *restrict s1;
    s1 = (Scalar_t *) malloc(sizeof(Scalar_t) * NUM_SPECIES
                                              * NUM_ESTEPS * NUM_MUSTEPS);

    for (s = 0; s < N_subcycles; s++) {

      if (AdiabaticChangeAlg == 3){
        AdiabaticChange_Operator_WENO3(f1,f,v_avg,node);
      } else {
        AdiabaticChange_Operator_Upwind(f1,f,v_avg,node);
      }

      for (species = 0; species < NUM_SPECIES; species++) {
        for (energy = 0; energy < NUM_ESTEPS; energy++) {
          for (mu = 0; mu < NUM_MUSTEPS; mu++) {

            idx = idx_spem(species,energy,mu);

            s1[idx] = f[idx] - dt_subcycle*f1[idx];

          }
        }
      }

      if (AdiabaticChangeAlg == 3){
        AdiabaticChange_Operator_WENO3(f1,s1,v_avg,node);
      } else {
        AdiabaticChange_Operator_Upwind(f1,s1,v_avg,node);
      }

      for (species = 0; species < NUM_SPECIES; species++) {
        for (energy = 0; energy < NUM_ESTEPS; energy++) {
          for (mu = 0; mu < NUM_MUSTEPS; mu++) {

            idx = idx_spem(species,energy,mu);

            s1[idx] = (3.0/4.0) * f[idx] + (1.0/4.0) * s1[idx]
                                         - (1.0/4.0) * dt_subcycle * f1[idx];

          }
        }
      }

      if (AdiabaticChangeAlg == 3){
        AdiabaticChange_Operator_WENO3(f1,s1,v_avg,node);
      } else {
        AdiabaticChange_Operator_Upwind(f1,s1,v_avg,node);
      }

      for (species = 0; species < NUM_SPECIES; species++) {
        for (energy = 0; energy < NUM_ESTEPS; energy++) {
          for (mu = 0; mu < NUM_MUSTEPS; mu++) {

            idx = idx_spem(species,energy,mu);

            f[idx] = (1.0/3.0) * f[idx] + (2.0/3.0) * s1[idx]
                                        - (2.0/3.0) * dt_subcycle * f1[idx];

          }
        }
      }
    } // subcycles

    free(s1);

  } else {

    for (s = 0; s < N_subcycles; s++) {

      AdiabaticChange_Operator_Upwind(f1,f,v_avg,node);

      for (species = 0; species < NUM_SPECIES; species++) {
        for (energy = 0; energy < NUM_ESTEPS; energy++) {
          for (mu = 0; mu < NUM_MUSTEPS; mu++) {

            idx = idx_spem(species,energy,mu);

            f[idx] = f[idx] - dt_subcycle * f1[idx];

          }
        }
      }

    } // subcycles
  }


  // Convert back to linear space and check for badness:

  for (species = 0; species < NUM_SPECIES; species++) {
    for (energy = 0; energy < NUM_ESTEPS; energy++) {
      for (mu = 0; mu < NUM_MUSTEPS; mu++) {

        idx = idx_frcsspem(face,row,col,shell,species,energy,mu);

        eParts[idx] = exp(f[idx_spem(species,energy,mu)]);

        // check for NaNs and Infs
        checkNaN(mpi_rank, face, row, col, shell,eParts[idx],
                 "Adiabatic Change");
        checkInf(mpi_rank, face, row, col, shell,eParts[idx],
                 "Adiabatic Change");

        // Check if distribution dropped below double minimum.
        // If so, set it to double minimum.
        if ( eParts[idx] < DBL_MIN ){
       //   warn(face, row, col, shell, species, energy, mu,
       //       "AdiabaticChange: distribution less than DBL_MIN", &eParts[idx]);
          eParts[idx] = DBL_MIN;
        }

      }
    }
  }

  // Free up temporary arrays.
  free(v_avg);
  free(vel);
  free(f1);
  free(f);

}
/*---------- END AdiabaticChange( ) --------------------------------*/
/*------------------------------------------------------------------*/

/*---------------------------------------------------------------*/
/*---------------------------------------------------------------*/
/*--*/    void                                               /*--*/
/*--*/    AdiabaticChange_Operator_Upwind(Scalar_t* f1,      /*--*/
/*--*/                                    Scalar_t* f,       /*--*/
/*--*/                                    Scalar_t* v_avg,   /*--*/
/*--*/                                    Node_t node)       /*--*/
/*--*/                                                       /*--*/
/*-------------------------------------------------------------- */
{

  Scalar_t upwind = 1.0;

  Scalar_t *restrict flux;
  Index_t species, energy, mu;
  Scalar_t cce;

  flux = (Scalar_t *) malloc(sizeof(Scalar_t) * NUM_SPECIES * (NUM_ESTEPS+1) * NUM_MUSTEPS);

  // First, compute "half-mesh" fluxes.

  for (species = 0; species < NUM_SPECIES; species++) {
    for (energy = 1; energy < NUM_ESTEPS; energy++) {
      for (mu = 0; mu < NUM_MUSTEPS; mu++) {

        cce = copysign(upwind,v_avg[idx_spep1m(species,energy,mu)]);
        flux[idx_spep1m(species,energy,mu)] = v_avg[idx_spep1m(species,energy,mu)]
           * 0.5 * ( (1.0 - cce) * f[idx_spem(species,energy  ,mu)]
                 +   (1.0 + cce) * f[idx_spem(species,energy-1,mu)]);

      }
    }
  }

  // Now update all internal points.

  for (species = 0; species < NUM_SPECIES; species++) {
    for (energy = 1; energy < NUM_ESTEPS-1; energy++) {
      for (mu = 0; mu < NUM_MUSTEPS; mu++) {

        f1[idx_spem(species,energy,mu)]=(flux[idx_spep1m(species,energy+1,mu)] -
                                         flux[idx_spep1m(species,energy  ,mu)])/dlnp;

      }
    }
  }

  // Boundary conditions:
  // For outflow, just use UW as usual.
  // For inflow, use Dirichlet of seed population for point "past the grid"

  for (species = 0; species < NUM_SPECIES; species++) {
    for (mu = 0; mu < NUM_MUSTEPS; mu++) {

      // Left boundary

      energy = 0;

      if (v_avg[idx_spep1m(species,energy+1,mu)] <= 0) {
        f1[idx_spem(species,energy,mu)] = v_avg[idx_spep1m(species,energy+1,mu)]
                                        *(f[idx_spem(species,energy+1,mu)] -
                                          f[idx_spem(species,energy  ,mu)])/dlnp;
      }
      else
      {
        f1[idx_spem(species,energy,mu)] = v_avg[idx_spep1m(species,energy+1,mu)]
                          *(f[idx_spem(species,energy,mu)] -
                    log(sepSeedFunction(egrid[idx_se(species, energy)], node.rmag)))/dlnp;
      }

      // Right boundary

      energy = NUM_ESTEPS-1;

      if (v_avg[idx_spep1m(species,energy,mu)] >= 0) {
        f1[idx_spem(species,energy,mu)] = v_avg[idx_spep1m(species,energy,mu)]
                                         *(f[idx_spem(species,energy  ,mu)] -
                                           f[idx_spem(species,energy-1,mu)])/dlnp;
      }
      else
      {
        f1[idx_spem(species,energy,mu)] = v_avg[idx_spep1m(species,energy,mu)]*
                    (log(sepSeedFunction(egrid[idx_se(species, energy)], node.rmag)) -
                                           f[idx_spem(species,energy,mu)])/dlnp;
      }
    }
  }

  free(flux);

} /*-------- END AdiabaticChange_Operator_Upwind()------------------*/
/*------------------------------------------------------------------*/

/*---------------------------------------------------------------*/
/*---------------------------------------------------------------*/
/*--*/    void                                               /*--*/
/*--*/    AdiabaticChange_Operator_WENO3(Scalar_t* f1,       /*--*/
/*--*/                                   Scalar_t* f,        /*--*/
/*--*/                                   Scalar_t* v_avg,    /*--*/
/*--*/                                   Node_t node )       /*--*/
/*--*/                                                       /*--*/
/*--  Does Weno3 and returns recompute vector                  --*/
/*-------------------------------------------------------------- */
{/*--------------------------------------------------------------*/

  Scalar_t *restrict LP, *restrict LN, *restrict alpha, *restrict flux;

  Scalar_t upwind = 1.0;

  Scalar_t weno_eps;
  Scalar_t D_C_CPt;
  Scalar_t D_C_MCt;
  Scalar_t D_M_Tt;
  Scalar_t D_CP_Tt;
  Scalar_t D_P_Tt;
  Scalar_t D_MC_Tt;
  Scalar_t vel_abs_max = 0.0;
  Scalar_t cce;

  Scalar_t p0m, p1m, p0p, p1p;
  Scalar_t B0m, B1m, B0p, B1p;
  Scalar_t w0m, w1m, w0p, w1p;
  Scalar_t wm_sum, wp_sum;
  Scalar_t OM0m, OM1m, OM0p, OM1p;
  Scalar_t um, up;
  Index_t species, energy, mu;

  D_C_CPt = 0.5;
  D_C_MCt = 0.5;
  D_M_Tt  = 1.0/3.0;
  D_CP_Tt = 2.0/3.0;
  D_P_Tt  = 1.0/3.0;
  D_MC_Tt = 2.0/3.0;
  weno_eps = 10.0*sqrt(DBL_MIN);

  flux  = (Scalar_t *) malloc(sizeof(Scalar_t) * NUM_SPECIES * (NUM_ESTEPS+1) * NUM_MUSTEPS);
  LP    = (Scalar_t *) malloc(sizeof(Scalar_t) * NUM_SPECIES * (NUM_ESTEPS) * NUM_MUSTEPS);
  LN    = (Scalar_t *) malloc(sizeof(Scalar_t) * NUM_SPECIES * (NUM_ESTEPS) * NUM_MUSTEPS);
  alpha = (Scalar_t *) malloc(sizeof(Scalar_t) * NUM_SPECIES * (NUM_ESTEPS) * NUM_MUSTEPS);

/* ****** Get alpha ****** */

/* NOTE: v_avg is on extended "half" mesh.
   This means it has a length in energy of NUM_ESTEPS+1
   Indicies:
   X    O    X    O    X    O    X    O    X    O    X
   0    0    1    1    2    2    3    3    NE-1 NE-1 NE
   (v_avg,flux on X), (alpha,LP/N,f,f1 on O)
*/

/* ***** Get inner points. */
  for (species = 0; species < NUM_SPECIES; species++) {
    for (energy = 2; energy < NUM_ESTEPS-2; energy++) {
      for (mu = 0; mu < NUM_MUSTEPS; mu++) {
        vel_abs_max =
            fabs(v_avg[idx_spep1m(species,energy  ,mu)]);
        if (fabs(v_avg[idx_spep1m(species,energy-1,mu)]) > vel_abs_max) vel_abs_max = fabs(v_avg[idx_spep1m(species,energy-1,mu)]);
        if (fabs(v_avg[idx_spep1m(species,energy-2,mu)]) > vel_abs_max) vel_abs_max = fabs(v_avg[idx_spep1m(species,energy-2,mu)]);
        if (fabs(v_avg[idx_spep1m(species,energy+1,mu)]) > vel_abs_max) vel_abs_max = fabs(v_avg[idx_spep1m(species,energy+1,mu)]);
        if (fabs(v_avg[idx_spep1m(species,energy+2,mu)]) > vel_abs_max) vel_abs_max = fabs(v_avg[idx_spep1m(species,energy+2,mu)]);
        if (fabs(v_avg[idx_spep1m(species,energy+3,mu)]) > vel_abs_max) vel_abs_max = fabs(v_avg[idx_spep1m(species,energy+3,mu)]);

        alpha[idx_spem(species,energy,mu)] = vel_abs_max;
      }
    }
  }

  /* ***** Edge cases. */

  for (species = 0; species < NUM_SPECIES; species++) {
    for (mu = 0; mu < NUM_MUSTEPS; mu++) {
      energy = 0;
      vel_abs_max =
          fabs(v_avg[idx_spep1m(species,energy  ,mu)]);
      if (fabs(v_avg[idx_spep1m(species,energy+1,mu)]) > vel_abs_max) vel_abs_max = fabs(v_avg[idx_spep1m(species,energy+1,mu)]);
      if (fabs(v_avg[idx_spep1m(species,energy+2,mu)]) > vel_abs_max) vel_abs_max = fabs(v_avg[idx_spep1m(species,energy+2,mu)]);
      if (fabs(v_avg[idx_spep1m(species,energy+3,mu)]) > vel_abs_max) vel_abs_max = fabs(v_avg[idx_spep1m(species,energy+3,mu)]);

      alpha[idx_spem(species,energy,mu)] = vel_abs_max;

      energy = 1;
      vel_abs_max =
          fabs(v_avg[idx_spep1m(species,energy  ,mu)]);
      if (fabs(v_avg[idx_spep1m(species,energy-1,mu)]) > vel_abs_max) vel_abs_max = fabs(v_avg[idx_spep1m(species,energy-1,mu)]);
      if (fabs(v_avg[idx_spep1m(species,energy+1,mu)]) > vel_abs_max) vel_abs_max = fabs(v_avg[idx_spep1m(species,energy+1,mu)]);
      if (fabs(v_avg[idx_spep1m(species,energy+2,mu)]) > vel_abs_max) vel_abs_max = fabs(v_avg[idx_spep1m(species,energy+2,mu)]);
      if (fabs(v_avg[idx_spep1m(species,energy+3,mu)]) > vel_abs_max) vel_abs_max = fabs(v_avg[idx_spep1m(species,energy+3,mu)]);

      alpha[idx_spem(species,energy,mu)] = vel_abs_max;

      energy = NUM_ESTEPS-2;
      vel_abs_max =
          fabs(v_avg[idx_spep1m(species,energy  ,mu)]);
      if (fabs(v_avg[idx_spep1m(species,energy+1,mu)]) > vel_abs_max) vel_abs_max = fabs(v_avg[idx_spep1m(species,energy+1,mu)]);
      if (fabs(v_avg[idx_spep1m(species,energy+2,mu)]) > vel_abs_max) vel_abs_max = fabs(v_avg[idx_spep1m(species,energy+2,mu)]);
      if (fabs(v_avg[idx_spep1m(species,energy-1,mu)]) > vel_abs_max) vel_abs_max = fabs(v_avg[idx_spep1m(species,energy-1,mu)]);
      if (fabs(v_avg[idx_spep1m(species,energy-2,mu)]) > vel_abs_max) vel_abs_max = fabs(v_avg[idx_spep1m(species,energy-2,mu)]);

      alpha[idx_spem(species,energy,mu)] = vel_abs_max;

      energy = NUM_ESTEPS-1;
      vel_abs_max =
          fabs(v_avg[idx_spep1m(species,energy  ,mu)]);
      if (fabs(v_avg[idx_spep1m(species,energy+1,mu)]) > vel_abs_max) vel_abs_max = fabs(v_avg[idx_spep1m(species,energy+1,mu)]);
      if (fabs(v_avg[idx_spep1m(species,energy-1,mu)]) > vel_abs_max) vel_abs_max = fabs(v_avg[idx_spep1m(species,energy-1,mu)]);
      if (fabs(v_avg[idx_spep1m(species,energy-2,mu)]) > vel_abs_max) vel_abs_max = fabs(v_avg[idx_spep1m(species,energy-2,mu)]);

      alpha[idx_spem(species,energy,mu)] = vel_abs_max;
    }
  }

  for (species = 0; species < NUM_SPECIES; species++) {
    for (energy = 0; energy < NUM_ESTEPS; energy++) {
      for (mu = 0; mu < NUM_MUSTEPS; mu++) {

        LP[idx_spem(species,energy,mu)] = 0.5 * f[idx_spem(species,energy,mu)] *
          v_avg[idx_spep1m(species,energy+1,mu)] - alpha[idx_spem(species,energy,mu)];

        LN[idx_spem(species,energy,mu)] = 0.5 * f[idx_spem(species,energy,mu)] *
          v_avg[idx_spep1m(species,energy,mu)] + alpha[idx_spem(species,energy,mu)];

      }
    }
  }

  for (species = 0; species < NUM_SPECIES; species++) {
    for (energy = 2; energy < NUM_ESTEPS-1; energy++) {
      for (mu = 0; mu < NUM_MUSTEPS; mu++) {

        p0m = (1.0 + D_C_MCt)*LN[idx_spem(species,energy-1,mu)] - D_C_MCt*LN[idx_spem(species,energy-2,mu)];
        p1m =         D_C_MCt *LN[idx_spem(species,energy-1,mu)] + D_C_CPt*LN[idx_spem(species,energy  ,mu)];
        p0p = (1.0 + D_C_CPt)*LP[idx_spem(species,energy  ,mu)] - D_C_CPt*LP[idx_spem(species,energy+1,mu)];
        p1p =         D_C_CPt *LP[idx_spem(species,energy  ,mu)] + D_C_MCt*LP[idx_spem(species,energy-1,mu)];

        B0m = 4.0*pow((D_C_MCt*(LN[idx_spem(species,energy-1,mu)] - LN[idx_spem(species,energy-2,mu)])),2);
        B1m = 4.0*pow((D_C_CPt*(LN[idx_spem(species,energy  ,mu)] - LN[idx_spem(species,energy-1,mu)])),2);
        B0p = 4.0*pow((D_C_CPt*(LP[idx_spem(species,energy+1,mu)] - LP[idx_spem(species,energy  ,mu)])),2);
        B1p = 4.0*pow((D_C_MCt*(LP[idx_spem(species,energy  ,mu)] - LP[idx_spem(species,energy-1,mu)])),2);


        w0m = D_P_Tt /pow((weno_eps + B0m),2);
        w1m = D_MC_Tt/pow((weno_eps + B1m),2);
        w0p = D_M_Tt /pow((weno_eps + B0p),2);
        w1p = D_CP_Tt/pow((weno_eps + B1p),2);

        wm_sum = w0m + w1m;
        wp_sum = w0p + w1p;

        OM0m = w0m/wm_sum;
        OM1m = w1m/wm_sum;
        OM0p = w0p/wp_sum;
        OM1p = w1p/wp_sum;

        um = OM0m*p0m + OM1m*p1m;
        up = OM0p*p0p + OM1p*p1p;

        flux[idx_spep1m(species,energy,mu)] = up + um;

      }
    }
  }

  // Use upwinding for fluxes near the boundary.

  for (species = 0; species < NUM_SPECIES; species++) {
    for (mu = 0; mu < NUM_MUSTEPS; mu++) {

      energy = 1;

      cce = copysign(upwind,v_avg[idx_spep1m(species,energy,mu)]);
      flux[idx_spep1m(species,energy,mu)] = v_avg[idx_spep1m(species,energy,mu)]
             * 0.5 * ( (1.0 - cce) * f[idx_spem(species,energy  ,mu)]
                   +   (1.0 + cce) * f[idx_spem(species,energy-1,mu)]);

      energy = NUM_ESTEPS-1;

      cce = copysign(upwind,v_avg[idx_spep1m(species,energy,mu)]);
      flux[idx_spep1m(species,energy,mu)] = v_avg[idx_spep1m(species,energy,mu)]
             * 0.5 * ( (1.0 - cce) * f[idx_spem(species,energy  ,mu)]
                   +   (1.0 + cce) * f[idx_spem(species,energy-1,mu)]);

    }
  }

  for (species = 0; species < NUM_SPECIES; species++) {
    for (energy = 1; energy < NUM_ESTEPS-1; energy++) {
      for (mu = 0; mu < NUM_MUSTEPS; mu++) {

        f1[idx_spem(species,energy,mu)]=(flux[idx_spep1m(species,energy+1,mu)] -
                                         flux[idx_spep1m(species,energy  ,mu)])/dlnp;

      }
    }
  }

  // Boundary conditions:
  // For outflow, just use UW as usual.
  // For inflow, use Dirichlet of seed population for point "past the grid"

  for (species = 0; species < NUM_SPECIES; species++) {
    for (mu = 0; mu < NUM_MUSTEPS; mu++) {

      // Left boundary

      energy = 0;

      if (v_avg[idx_spep1m(species,energy+1,mu)] <= 0) {
        f1[idx_spem(species,energy,mu)] = v_avg[idx_spep1m(species,energy+1,mu)]
                                        *(f[idx_spem(species,energy+1,mu)] -
                                          f[idx_spem(species,energy  ,mu)])/dlnp;
      }
      else
      {
        f1[idx_spem(species,energy,mu)] = v_avg[idx_spep1m(species,energy+1,mu)]
                          *(f[idx_spem(species,energy,mu)] -
                    log(sepSeedFunction(egrid[idx_se(species, energy)], node.rmag)))/dlnp;
      }

      // Right boundary

      energy = NUM_ESTEPS-1;

      if (v_avg[idx_spep1m(species,energy,mu)] >= 0) {
        f1[idx_spem(species,energy,mu)] = v_avg[idx_spep1m(species,energy,mu)]
                                         *(f[idx_spem(species,energy  ,mu)] -
                                           f[idx_spem(species,energy-1,mu)])/dlnp;
      }
      else
      {
        f1[idx_spem(species,energy,mu)] = v_avg[idx_spep1m(species,energy,mu)]*
                    (log(sepSeedFunction(egrid[idx_se(species, energy)], node.rmag)) -
                                           f[idx_spem(species,energy,mu)])/dlnp;
      }
    }
  }

  free(flux);
  free(LP);
  free(LN);
  free(alpha);

} /*-------- END AdiabaticChange_Operator_WENO3( ) -----------------*/
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/*--*/             void                                         /*--*/
/*--*/     ShockSolution( Index_t face,                         /*--*/
/*--*/                    Index_t row,                          /*--*/
/*--*/                    Index_t col,                          /*--*/
/*--*/                    Index_t shell,                        /*--*/
/*--*/                    Index_t species,                      /*--*/
/*--*/                    Scalar_t dt)                          /*--*/
/*--                                                              --*/
/*-- Evaluate the adiabatic change                                --*/
/*------------------------------------------------------------------*/
{/*-----------------------------------------------------------------*/

  Node_t node;
  Scalar_t r0, r1, b0, b1, u0, u1, n0, n1;
  Scalar_t cr, gamma, beta;
  Scalar_t b_dn2xb_up2, cos2_tBN, sin2_tBN, thetaBN, denominator;
  Scalar_t vInj, pInj, eInj;
  Scalar_t pMin, pMax, pProj, pProjDist, pProjStepGrid;
  Scalar_t p, p1, p2, f1, f2;
  Index_t mu, energy, eMinStep, fInjIndex, pProjStep;
  Scalar_t fInj, p_fInjIndex;
  Scalar_t oneOverGamma, v, lambda0, lambda1;
  Scalar_t kappa0, kappa1, kappaPar0, kappaPar1, kappaPerp0, kappaPerp1;
  Scalar_t ionGyroFreq0, ionGyroFreq1, ionGyroRadius0, ionGyroRadius1;
  Scalar_t dx, deltaU;
  Scalar_t seedFunction;

  node = grid[idx_frcs(face,row,col,shell)];


  // the compiler throws up a complaint about these maybe being uninitialized, so I'm setting a value here to keep
  // it quiet.  If they're used, they are initialized, I've gone through the logic.
  p_fInjIndex = 0;
  fInj = 0.0;
  fInjIndex = 0;

  // paramaters necessary for shock solution

  // there is a slight funny business w/ regards to the up and downstream positions:
  // due to everything being measured in the moving frame, the position of the downstream
  // node ends up being ahead of the upstream node.  The only place this is used is for the
  // seed function and the spatial difference is so small it probably doesn't matter.  But,
  // just to be consistent, we're swapping their positions.
  r0 = node.rmag;
  r1 = sqrt(node.rOlder.x*node.rOlder.x + node.rOlder.y*node.rOlder.y + node.rOlder.z*node.rOlder.z);

  n1 = node.mhdDensity;
  n0 = node.mhdDensityOld;
  b1 = node.mhdBmag;
  b0 = node.mhdBmagOld;
  u1 = node.mhdVmag;
  u0 = sqrt(node.mhdVsphOld.r*node.mhdVsphOld.r + node.mhdVsphOld.theta*node.mhdVsphOld.theta + node.mhdVsphOld.phi*node.mhdVsphOld.phi);

  // compression ratio and check against maximum
  cr = n1 / n0;

  if (cr > 4.0)
    cr = 4.0;

  gamma = 3.0 * cr / (cr - 1.0);

  if (cr > 1.3) {

    // deltaU
    deltaU = sqrt(pow(node.r.x - node.rOlder.x,2.0) + pow(node.r.y - node.rOlder.y,2.0) + pow(node.r.z - node.rOlder.z,2.0)) * config.rScale * log(cr);

    // determining the shock angle and injection energy
    b_dn2xb_up2 = b1 * b1 / (b0 * b0 + VERYSMALL);

    if ( fabs(b_dn2xb_up2 - 1.0) > (cr * cr - 1.0) ) {

      cos2_tBN = 0.0;
      sin2_tBN = 1.0;
      thetaBN = PI / 2.0;
      vInj = u0;

    } else {

      thetaBN = asin( sqrt(fabs(b_dn2xb_up2 - 1.0) / (cr * cr - 1.0)) );
      sin2_tBN = sin(thetaBN) * sin(thetaBN);
      cos2_tBN = cos(thetaBN) * cos(thetaBN);
      denominator = config.kperxkpar * sin2_tBN + cos2_tBN;
      vInj = u0 * sqrt( 1.0 + ((1.0 - config.kperxkpar) * (1.0 - config.kperxkpar) * sin2_tBN * cos2_tBN /
                               (denominator * denominator)) );

    }

    // injection momentum, energy, and distribution
    pInj = vInj / (sqrt(1.0 - vInj * vInj) + VERYSMALL);
    eInj = sqrt( 1.0 + pInj * pInj ) - 1.0;

    // making sure there is a lower bound on the injection energy
    if (eInj < (config.minInjectionEnergy * MEV / (MP * C * C))) {

      eInj = config.minInjectionEnergy * MEV / (MP * C * C);
      pInj = sqrt((eInj + 1.0) * (eInj + 1.0) - 1.0);

    }

    // storing the seed function, it is called many times later
    seedFunction = sepSeedFunction(eInj, r0) * config.shockInjectionFactor;

    // determining the lowest step on the grid to be accelerated and the injection
    pMin = exp(lnpmin);
    pMax = exp(lnpmin + (NUM_ESTEPS - 1) * dlnp);

    if (pInj < pMin) {

      fInj = seedFunction;
      eMinStep = 0;

    } else if (pInj >= pMax) {

      eMinStep = NUM_ESTEPS;

    } else {

      fInjIndex = floor( log(pInj / pMin) / dlnp + 0.5);
      p_fInjIndex = exp(lnpmin + fInjIndex * dlnp);

      if (fInjIndex >= (NUM_ESTEPS - 1))
        eMinStep = NUM_ESTEPS;
      else
        eMinStep = fInjIndex + 1;

    }

    ionGyroFreq1 = fabs(config.charge[species]) * OM * b1 / config.mass[species];
    ionGyroFreq0 = fabs(config.charge[species]) * OM * b0 / config.mass[species];

    // looping over the energy bins which get accelerated
    for (energy = eMinStep; energy < NUM_ESTEPS; energy++) {

      // determine the scale length
      v = vgrid[energy];
      p = exp(lnpmin + energy * dlnp);

      lambda1 = meanFreePath(species, energy, r1 * config.rScale);
      lambda0 = meanFreePath(species, energy, r0 * config.rScale);

      kappaPar1 = v * lambda1 / 3.0;
      kappaPar0 = v * lambda0 / 3.0;

      oneOverGamma = v / p;

      ionGyroRadius1 = v / ( oneOverGamma * ionGyroFreq1 + VERYSMALL);
      ionGyroRadius0 = v / ( oneOverGamma * ionGyroFreq0 + VERYSMALL);

      kappaPerp1 = kappaPar1 / ( 1.0 + (lambda1 * lambda1) / (ionGyroRadius1 * ionGyroRadius1 + VERYSMALL) );
      kappaPerp0 = kappaPar0 / ( 1.0 + (lambda0 * lambda0) / (ionGyroRadius0 * ionGyroRadius0 + VERYSMALL) );

      kappa1 = kappaPar1 * cos2_tBN + kappaPerp1 * sin2_tBN;
      kappa0 = kappaPar0 * cos2_tBN + kappaPerp0 * sin2_tBN;

      dx = kappa1 / (u1 + VERYSMALL) + kappa0 / (u0 + VERYSMALL);

      // calculated the backward projected momentum
      pProj = p * exp(-1.0 * dt * deltaU / (3.0 * dx + VERYSMALL));
      pProjStep = floor( log(pProj/pMin) / dlnp + 0.5 );

      if (pProjStep < energy) {

        for (mu = 0; mu < NUM_MUSTEPS; mu++) {

          // determining the injection distribution if the injection energy is on the grid
          // i'm hardcoding in the logic but it should probably be made a function
          if (pInj >= pMin) {

            if (pInj > p_fInjIndex) {

              p1 = pgrid[fInjIndex];
              f1 = eParts[idx_frcsspem(face,row,col,shell,species,fInjIndex,mu)];
              p2 = pgrid[fInjIndex+1];
              f2 = eParts[idx_frcsspem(face,row,col,shell,species,fInjIndex+1,mu)];

              if ((f1 == 0.0) || (f2 == 0.0))
                fInj = linInterp(f1, f2, pInj, p1, p2);
              else {

                beta = log(f2 / f1) / dlnp;

                if (fabs(beta) > 1.0)
                  fInj = f1 * pow(pInj / p_fInjIndex, beta);
                else
                  fInj = linInterp(f1, f2, pInj, p1, p2);

              }

            } else if (pInj < p_fInjIndex) {

              p1 = pgrid[fInjIndex-1];
              f1 = eParts[idx_frcsspem(face,row,col,shell,species,fInjIndex-1,mu)];
              p2 = pgrid[fInjIndex];
              f2 = eParts[idx_frcsspem(face,row,col,shell,species,fInjIndex,mu)];

              if ((f1 == 0.0) || (f2 == 0.0))
                fInj = linInterp(f1, f2, pInj, p1, p2);
              else {

                beta = log(f2 / f1) / dlnp;

                if (fabs(beta) > 1.0)
                  fInj = f1 * pow(p_fInjIndex / pInj, beta);
                else
                  fInj = linInterp(f1, f2, pInj, p1, p2);

              }

            } else
              fInj = eParts[idx_frcsspem(face,row,col,shell,species,fInjIndex,mu)];

            if (fInj < seedFunction)
              fInj = seedFunction;

          }

          //if (mpi_rank == 2)
          //  printf("e:%i\tcr:%0.3f\tpMin:%0.3e\tpInj:%0.3e\tpProj:%0.3e\tpProjStep:%i\n", energy, cr, pMin, pInj, pProj, pProjStep);

          // three cases: pProj < pInj; pInj < pProj < pMin; pMin < pProj
          if (pProj <= pInj) {

            shockDist[idx_frcsspem(face,row,col,shell,species,energy,mu)] = fInj * pow(p / pInj, -1.0 * gamma);
            //eParts[idx_frcsspem(face,row,col,shell,species,energy,mu)] = fInj * pow(p / pInj, -1.0 * gamma);

          } else if (pProj < pMin) {

            p1 = pInj;
            f1 = fInj;
            p2 = pMin;
            f2 = eParts[idx_frcsspem(face,row,col,shell,species,0,mu)];

            if ((f1 == 0.0) || (f2 == 0.0))
              shockDist[idx_frcsspem(face,row,col,shell,species,energy,mu)] = linInterp(f1, f2, pProj, p1, p2);
              //eParts[idx_frcsspem(face,row,col,shell,species,energy,mu)] = linInterp(f1, f2, pProj, p1, p2);

            else {

              beta = log(f2 / f1) / log(p2 / p1);

              if (fabs(beta) > 1.0)
                shockDist[idx_frcsspem(face,row,col,shell,species,energy,mu)] = fInj * pow(pProj / pInj, beta) * pow(p / pProj, -1.0 * gamma);
                //eParts[idx_frcsspem(face,row,col,shell,species,energy,mu)] = fInj * pow(pProj / pInj, beta) * pow(p / pProj, -1.0 * gamma);
              else
                shockDist[idx_frcsspem(face,row,col,shell,species,energy,mu)] = linInterp(f1, f2, pProj, p1, p2);
                //eParts[idx_frcsspem(face,row,col,shell,species,energy,mu)] = linInterp(f1, f2, pProj, p1, p2);

            }

          } else {

            if (pProjStep <= 0) {

              pProj = pMin;
              pProjDist = eParts[idx_frcsspem(face,row,col,shell,species,0,mu)];

              if (pProjDist < (sepSeedFunction(egrid[idx_se(species,0)], r0) * config.shockInjectionFactor))
                pProjDist = sepSeedFunction(egrid[idx_se(species,0)], r0) * config.shockInjectionFactor;

              shockDist[idx_frcsspem(face,row,col,shell,species,energy,mu)] = pProjDist * pow(p / pProj, -1.0 * gamma);
              //eParts[idx_frcsspem(face,row,col,shell,species,energy,mu)] = pProjDist * pow(p / pProj, -1.0 * gamma);

            } else if (pProjStep >= (NUM_ESTEPS - 1)) {

              pProj = pMax;
              pProjDist = eParts[idx_frcsspem(face,row,col,shell,species,NUM_ESTEPS - 1,mu)];

              if (pProjDist < (sepSeedFunction(egrid[idx_se(species,NUM_ESTEPS - 1)], r0) * config.shockInjectionFactor))
                pProjDist = sepSeedFunction(egrid[idx_se(species,NUM_ESTEPS - 1)], r0) * config.shockInjectionFactor;

              shockDist[idx_frcsspem(face,row,col,shell,species,energy,mu)] = pProjDist * pow(p / pProj, -1.0 * gamma);
              //eParts[idx_frcsspem(face,row,col,shell,species,energy,mu)] = pProjDist * pow(p / pProj, -1.0 * gamma);

            } else {

              pProjStepGrid = pgrid[pProjStep];

              if (pProj > pProjStepGrid) {

                p1 = pgrid[pProjStep];
                f1 = eParts[idx_frcsspem(face,row,col,shell,species,pProjStep,mu)];
                p2 = pgrid[pProjStep+1];
                f2 = eParts[idx_frcsspem(face,row,col,shell,species,pProjStep+1,mu)];

                if ((f1 == 0.0) || (f2 == 0.0))
                  pProjDist = linInterp(f1, f2, pProj, p1, p2);
                else {

                  beta = log(f2 / f1) / dlnp;

                  if (fabs(beta) > 1.0)
                    pProjDist = f1 * pow(pProj / pProjStepGrid, beta);
                  else
                    pProjDist = linInterp(f1, f2, pProj, p1, p2);

                }

              } else if (pProj < pProjStepGrid) {

                p1 = pgrid[pProjStep-1];
                f1 = eParts[idx_frcsspem(face,row,col,shell,species,pProjStep - 1,mu)];
                p2 = pgrid[pProjStep];
                f2 = eParts[idx_frcsspem(face,row,col,shell,species,pProjStep,mu)];

                if ((f1 == 0.0) || (f2 == 0.0))
                  pProjDist = linInterp(f1, f2, pProj, p1, p2);
                else {

                  beta = log(f2 / f1) / dlnp;

                  if (fabs(beta) > 1.0)
                    pProjDist = f1 * pow(pProjStepGrid / pProj, beta);
                  else
                    pProjDist = linInterp(f1, f2, pProj, p1, p2);

                }

              } else
                pProjDist = eParts[idx_frcsspem(face,row,col,shell,species,pProjStep,mu)];

              if (pProjDist < (sepSeedFunction(egrid[idx_se(species,pProjStep)], r0) * config.shockInjectionFactor))
                pProjDist = sepSeedFunction(egrid[idx_se(species,pProjStep)], r0) * config.shockInjectionFactor;

              shockDist[idx_frcsspem(face,row,col,shell,species,energy,mu)] = pProjDist * pow(p / pProj, -1.0 * gamma);
              //eParts[idx_frcsspem(face,row,col,shell,species,energy,mu)] = pProjDist * pow(p / pProj, -1.0 * gamma);

            }

          }

          if (shockDist[idx_frcsspem(face,row,col,shell,species,energy,mu)] < eParts[idx_frcsspem(face,row,col,shell,species,energy,mu)])
            shockDist[idx_frcsspem(face,row,col,shell,species,energy,mu)] = eParts[idx_frcsspem(face,row,col,shell,species,energy,mu)];

        }

      }

    }

  }

}
/*---------- END ShockSolver( ) ------------------------------------*/
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/*--*/     void                                                 /*--*/
/*--*/     AdiabaticFocusing_old_mod( Index_t face,             /*--*/
/*--*/                                Index_t row,              /*--*/
/*--*/                                Index_t col,              /*--*/
/*--*/                                Index_t shell,            /*--*/
/*--*/                                Scalar_t dt )             /*--*/
/*--                                                              --*/
/*-- Calculate adiabatic focusing along a stream                  --*/
/*--  THIS SHOULD BE UPDATED TO USE ORIGINAL CHANGE ALG!          --*/
/*------------------------------------------------------------------*/
{/*-----------------------------------------------------------------*/

  Node_t node;

  Index_t species, energy, mu;

  Scalar_t dmu, mum, mup, muval;
  Scalar_t Af, Bf, Cf;
  Scalar_t beta, delta, Xi;
  Scalar_t dlnBds;
  Scalar_t* delEP;

  delEP = (Scalar_t*)malloc(sizeof(Scalar_t)*NUM_MUSTEPS);

  node = grid[idx_frcs(face,row,col,shell)];
//
// RMC: NOTE! These values (ds,Bmag+/-) are set
//            in update_stream_from_shell()+updateStreamValues()!
//
  if ((node.mhdBmagPlus == 0.0) || (node.mhdBmagMinus == 0.0))
    dlnBds = 0.0;
  else
    dlnBds = (log(node.mhdBmagPlus) - log(node.mhdBmagMinus)) / (2.0 * node.ds + DBL_MIN);

  Cf = (2.0 * node.mhdDlnN - 3.0 * node.mhdDlnB) / config.numEpSteps;

  for (species = 0; species < NUM_SPECIES; species++)
  {

    for (energy = 0; energy < NUM_ESTEPS; energy++)
    {

      Af = -vgrid[energy] * dlnBds * dt;
      Bf = (2.0 / vgrid[energy]) * node.mhdDuPar / (1.0 * config.numEpSteps);

      for (mu = 0; mu < NUM_MUSTEPS; mu++)
        delEP[mu] = 0.0;

      for (mu = 0; mu < NUM_MUSTEPS; mu++)
      {

        muval = mugrid[mu];

//        mum = -1.0 + mu * dmu;
//        mup = -1.0 + (mu + 1.0) * dmu;

//        beta = (Af - Bf) * (1.0 / (6.0 * dmu)) * (3.0 * dmu - mup*mup*mup + mum*mum*mum) +
//        Cf * (1.0 / (8.0 * dmu)) * (2.0 * (mup*mup - mum*mum) - mup*mup*mup*mup + mum*mum*mum*mum);

        beta = 0.5*(1.0-muval*muval)*((Af - Bf) + muval*Cf);

        delta = beta / dmu;

        if ( (delta > 0.0) && (mu < (NUM_MUSTEPS - 1)) )
        {
          delEP[mu]     -= delta * eParts[idx_frcsspem(face,row,col,shell,species,energy,mu)];
          delEP[mu + 1] += delta * eParts[idx_frcsspem(face,row,col,shell,species,energy,mu)];
        }

        if ( (delta < 0.0) && (mu > 0) )
        {
          delEP[mu]     += delta * eParts[idx_frcsspem(face,row,col,shell,species,energy,mu)];
          delEP[mu - 1] -= delta * eParts[idx_frcsspem(face,row,col,shell,species,energy,mu)];
        }

      }

      for (mu = 0; mu < NUM_MUSTEPS; mu++)
      {

        eParts[idx_frcsspem(face,row,col,shell,species,energy,mu)] += delEP[mu];

        // check for NaNs and Infs
        checkNaN(mpi_rank, face, row, col, shell,
                 eParts[idx_frcsspem(face,row,col,shell,species,energy,mu)],
                 "Adiabatic Focusing");
        checkInf(mpi_rank, face, row, col, shell,
                 eParts[idx_frcsspem(face,row,col,shell,species,energy,mu)],
                 "Adiabatic Focusing");

        // Check if distribution dropped below double minimum.
        // If so, set it to double minimum.
        // Note that this may make this scheme not conservative
        // but only by a tiny amount so should be OK.
        if ( eParts[idx_frcsspem(face,row,col,shell,species,energy,mu)] < DBL_MIN ){
          eParts[idx_frcsspem(face,row,col,shell,species,energy,mu)] = DBL_MIN;
        }
      }
    }
  }

  free(delEP);

}/*-------- END AdiabaticFocusing_old() --------------------------*/
/*------------------------------------------------------------------*/

/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/*--*/     void                                                 /*--*/
/*--*/     AdiabaticFocusing( Index_t face,                     /*--*/
/*--*/                        Index_t row,                      /*--*/
/*--*/                        Index_t col,                      /*--*/
/*--*/                        Index_t shell,                    /*--*/
/*--*/                        Scalar_t dt )                     /*--*/
/*--                                                              --*/
/*-- Calculate adiabatic focusing along a stream                  --*/
/*--                                                              --*/
/*------------------------------------------------------------------*/
{/*-----------------------------------------------------------------*/

// Want to solve f^n+1 = f^n - dt * vmu * df/dmu
// Here, all d/dt MHD derivative terms have not yet been divided by dt.
// The "dt" passed into this function is the EP subcyled time step.
// Therefore, we calcualate the large DT (the full EPREM time-step) to
// compute the MHD derivative constants.

  Node_t node;

  Index_t N_subcycles                      = 1;
  Scalar_t safety_factor                   = 0.9;
  Scalar_t vel_abs_max                     = 0.0;
  Scalar_t half                            = 0.5;
  Scalar_t one                             = 1.0;
  Scalar_t two                             = 2.0;
  Scalar_t three                           = 3.0;

  Index_t s, species, energy, mu, idx;

  Scalar_t *restrict vel, *restrict f, *restrict f1;
  Scalar_t *restrict v_avg;

  Scalar_t dt_full, dt_stable, dt_subcycle;
  Scalar_t DlnBDt, DlnNDt, DuParDt, dlnBds;
  Scalar_t a, b, muval;

  // Allocate space for the dist function, effective advection velocity, and fluxes.

  f    = (Scalar_t *) malloc(sizeof(Scalar_t) * NUM_SPECIES * NUM_ESTEPS * NUM_MUSTEPS);
  f1   = (Scalar_t *) malloc(sizeof(Scalar_t) * NUM_SPECIES * NUM_ESTEPS * NUM_MUSTEPS);
  vel  = (Scalar_t *) malloc(sizeof(Scalar_t) * NUM_SPECIES * NUM_ESTEPS * NUM_MUSTEPS);
  v_avg  = (Scalar_t *) malloc(sizeof(Scalar_t) * NUM_SPECIES * NUM_ESTEPS * (NUM_MUSTEPS+1));

  // Get the current node.  This contains the MHD differences
  // after the node has been moved (e.g. Delta-MHD = MHD^n+1-MHD^n)

  node = grid[idx_frcs(face,row,col,shell)];

  // Get the full timestep value and use it to compute MHD derivative terms.
  dt_full = dt*config.numEpSteps;
  DlnBDt  = node.mhdDlnB/dt_full;
  DlnNDt  = node.mhdDlnN/dt_full;
  DuParDt = node.mhdDuPar/dt_full;

  // Set spatial derivative of b-hat*Grad[ln(B)]:
  // NOTE! These values (ds,Bmag+/-) are set
  //       in update_stream_from_shell() + updateStreamValues()!
  // The use of DBL_MIN should be unnessesary as ds should never be zero except during
  // node seeding (in which case this routine should not be being called!).
  // Leaving it in for now...
  // ALSO, this is a central difference..   maybe it should be upwinded in the direction of B?
  if ((node.mhdBmagPlus == 0.0) || (node.mhdBmagMinus == 0.0))
    dlnBds = 0.0;
  else
    dlnBds = (log(node.mhdBmagPlus) - log(node.mhdBmagMinus)) / (two * node.ds + DBL_MIN);

  // Compute velocites.
  // Copy the distribution function into f
  // in case we want to later test integrating in ln space.

  b = two * DlnNDt - three * DlnBDt;

  for (species = 0; species < NUM_SPECIES; species++) {
    for (energy = 0; energy < NUM_ESTEPS; energy++) {
      for (mu = 0; mu < NUM_MUSTEPS; mu++) {

        idx = idx_spem(species,energy,mu);

        f[idx] = eParts[idx_frcsspem(face,row,col,shell,species,energy,mu)];

        muval = mugrid[mu];

        a = -vgrid[energy] * dlnBds - (two / vgrid[energy]) * DuParDt;

        vel[idx] = half*(one-muval*muval)*(a + muval*b);

      }
    }
  }

  // Generate extended v_avg array.
  // Eventually, and alternative vgrid array should be made to
  // match the "half mesh" directly so no averaging is needed.
  // Find max |v| for CFL condition.

  for (species = 0; species < NUM_SPECIES; species++) {
    for (energy = 0; energy < NUM_ESTEPS; energy++) {
      for (mu = 0; mu < NUM_MUSTEPS-1; mu++) {

        v_avg[idx_spemp1(species,energy,mu+1)] = half * (
                          vel[idx_spem(species,energy,mu+1)]
                        + vel[idx_spem(species,energy,mu  )]);

        if (fabs(v_avg[idx_spemp1(species,energy,mu+1)]) > vel_abs_max){
          vel_abs_max = fabs(v_avg[idx_spemp1(species,energy,mu+1)]);
        }
      }
    }
  }

  // Make extended v_avg boundary values (just copy for now).
  // Since copied, does not effect CFL.

  for (species = 0; species < NUM_SPECIES; species++) {
    for (energy = 0; energy < NUM_ESTEPS; energy++) {
      mu = 0;
      v_avg[idx_spemp1(species,energy,mu)] = v_avg[idx_spemp1(species,energy,mu+1)];
      mu = NUM_MUSTEPS;
      v_avg[idx_spemp1(species,energy,mu)] = v_avg[idx_spemp1(species,energy,mu-1)];
    }
  }

  // Find the stable time-step based on the velocity and dmu.

  dt_stable = safety_factor * dmu/vel_abs_max;

  // Now that we have the stable timestep, we need to be careful
  // about how many subcycles to use.
  // This routine is being called within the
  // EpSubCycle loop, so the dt we need to go is dt_full/config.numEpSteps,
  // or, simply the "dt" sent into the routine.

  N_subcycles = (Index_t) ceil(dt/dt_stable);

  // Keep track of the maximum number of subcycles
  if (N_subcycles > maxsubcycles_focusing){
    maxsubcycles_focusing = N_subcycles;
  }

  // Set the dt_subcycle so that it exactly reaches dt.

  dt_subcycle = dt/N_subcycles;

  // Now compute the sub-cycled upwind advance.
  if (AdiabaticFocusAlg == 2 || AdiabaticFocusAlg == 3){

    Scalar_t *restrict s1;
    s1 = (Scalar_t *) malloc(sizeof(Scalar_t) * NUM_SPECIES
                                              * NUM_ESTEPS * NUM_MUSTEPS);

    for (s = 0; s < N_subcycles; s++) {

      if (AdiabaticFocusAlg == 3){
        AdiabaticFocusing_Operator_WENO3(f1,f,v_avg);
      } else {
        AdiabaticFocusing_Operator_Upwind(f1,f,v_avg);
      }

      for (species = 0; species < NUM_SPECIES; species++) {
        for (energy = 0; energy < NUM_ESTEPS; energy++) {
          for (mu = 0; mu < NUM_MUSTEPS; mu++) {

            idx = idx_spem(species,energy,mu);

            s1[idx] = f[idx] - dt_subcycle*f1[idx];

          }
        }
      }

      if (AdiabaticFocusAlg == 3){
        AdiabaticFocusing_Operator_WENO3(f1,s1,v_avg);
      } else {
        AdiabaticFocusing_Operator_Upwind(f1,s1,v_avg);
      }

      for (species = 0; species < NUM_SPECIES; species++) {
        for (energy = 0; energy < NUM_ESTEPS; energy++) {
          for (mu = 0; mu < NUM_MUSTEPS; mu++) {

            idx = idx_spem(species,energy,mu);

            s1[idx] = (3.0/4.0) * f[idx] + (1.0/4.0) * s1[idx]
                                         - (1.0/4.0) * dt_subcycle * f1[idx];

          }
        }
      }

      if (AdiabaticFocusAlg == 3){
        AdiabaticFocusing_Operator_WENO3(f1,s1,v_avg);
      } else {
        AdiabaticFocusing_Operator_Upwind(f1,s1,v_avg);
      }

      for (species = 0; species < NUM_SPECIES; species++) {
        for (energy = 0; energy < NUM_ESTEPS; energy++) {
          for (mu = 0; mu < NUM_MUSTEPS; mu++) {

            idx = idx_spem(species,energy,mu);

            f[idx] = (1.0/3.0) * f[idx] + (2.0/3.0) * s1[idx]
                                        - (2.0/3.0) * dt_subcycle * f1[idx];

          }
        }
      }
    } // subcycles

    free(s1);

  } else {

    for (s = 0; s < N_subcycles; s++) {

      AdiabaticFocusing_Operator_Upwind(f1,f,v_avg);

      for (species = 0; species < NUM_SPECIES; species++) {
        for (energy = 0; energy < NUM_ESTEPS; energy++) {
          for (mu = 0; mu < NUM_MUSTEPS; mu++) {

            idx = idx_spem(species,energy,mu);

            f[idx] = f[idx] - dt_subcycle * f1[idx];

          }
        }
      }

    } // subcycles
  }

  // Copy back new distribution and check for badness:

  for (species = 0; species < NUM_SPECIES; species++) {
    for (energy = 0; energy < NUM_ESTEPS; energy++) {
      for (mu = 0; mu < NUM_MUSTEPS; mu++) {

        idx = idx_frcsspem(face,row,col,shell,species,energy,mu);

        eParts[idx] = f[idx_spem(species,energy,mu)];

        // check for NaNs and Infs
        checkNaN(mpi_rank, face, row, col, shell,eParts[idx],
                 "Adiabatic Focusing");
        checkInf(mpi_rank, face, row, col, shell,eParts[idx],
                 "Adiabatic Focusing");

        // Check if distribution dropped below double minimum.
        // If so, set it to double minimum.
        // Note that this may make this scheme not conservative
        // but only by a tiny amount so should be OK.
        if ( eParts[idx] < DBL_MIN ){
          eParts[idx] = DBL_MIN;
        }

      }
    }
  }

  // Free up temporary arrays.
  free(v_avg);
  free(vel);
  free(f1);
  free(f);

}/*-------- END AdiabaticFocusing() --------------------------------*/
/*------------------------------------------------------------------*/

/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/*--*/    void                                                  /*--*/
/*--*/    AdiabaticFocusing_Operator_Upwind(Scalar_t* f1,       /*--*/
/*--*/                                      Scalar_t* f,        /*--*/
/*--*/                                      Scalar_t* v_avg)    /*--*/
/*--*/                                                          /*--*/
/*----------------------------------------------------------------- */
{

  Scalar_t upwind = 1.0;

  Scalar_t *restrict flux;
  Index_t species, energy, mu;
  Scalar_t cce;
  Scalar_t dmu_i;

  dmu_i = 1.0/dmu;

  flux = (Scalar_t *) malloc(sizeof(Scalar_t) * NUM_SPECIES * NUM_ESTEPS * (NUM_MUSTEPS+1));

  // First, compute "half-mesh" inner fluxes.

  for (species = 0; species < NUM_SPECIES; species++) {
    for (energy = 0; energy < NUM_ESTEPS; energy++) {
      for (mu = 1; mu < NUM_MUSTEPS; mu++) {

      cce = copysign(upwind,v_avg[idx_spemp1(species,energy,mu)]);
      flux[idx_spemp1(species,energy,mu)] = v_avg[idx_spemp1(species,energy,mu)]
             * 0.5 * ( (1.0 - cce) * f[idx_spem(species,energy,mu  )]
                   +   (1.0 + cce) * f[idx_spem(species,energy,mu-1)]);

      }
    }
  }

  // Now update all internal points.

  for (species = 0; species < NUM_SPECIES; species++) {
    for (energy = 0; energy < NUM_ESTEPS; energy++) {
      for (mu = 1; mu < NUM_MUSTEPS-1; mu++) {

        f1[idx_spem(species,energy,mu)]=(flux[idx_spemp1(species,energy,mu+1)] -
                                         flux[idx_spemp1(species,energy,mu  )])*dmu_i;

      }
    }
  }

  // Boundary conditions:
  // Solid wall
  // For inflow, allow stuff to enter boundary cell, but not leave it.
  // For outflow, allow stuff to leave boundary, nothing enters.

  for (species = 0; species < NUM_SPECIES; species++) {
    for (energy = 0; energy < NUM_ESTEPS; energy++) {

      // Left boundary

      mu = 0;

      if (v_avg[idx_spemp1(species,energy,mu+1)] <= 0) {
        f1[idx_spem(species,energy,mu)] = v_avg[idx_spemp1(species,energy,mu+1)]
                                             *f[idx_spem(species,energy,mu+1)]*dmu_i;
      }
      else
      {
        f1[idx_spem(species,energy,mu)] = v_avg[idx_spemp1(species,energy,mu+1)]
                                             *f[idx_spem(species,energy,mu)]*dmu_i;
      }

      // Right boundary

      mu = NUM_MUSTEPS-1;

      if (v_avg[idx_spemp1(species,energy,mu)] >= 0) {
        f1[idx_spem(species,energy,mu)] = v_avg[idx_spemp1(species,energy,mu)]
                                           *(-f[idx_spem(species,energy,mu-1)])*dmu_i;
      }
      else
      {
        f1[idx_spem(species,energy,mu)] = v_avg[idx_spemp1(species,energy,mu)]
                                           *(-f[idx_spem(species,energy,mu)])*dmu_i;
      }
    }
  }

  free(flux);

}
/*-------- END AdiabaticFocusing_Operator_Upwind()------------------*/
/*------------------------------------------------------------------*/

/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/*--*/    void                                                  /*--*/
/*--*/    AdiabaticFocusing_Operator_WENO3(Scalar_t* f1,        /*--*/
/*--*/                                     Scalar_t* f,         /*--*/
/*--*/                                     Scalar_t* v_avg)     /*--*/
/*--*/                                                          /*--*/
/*--  Does Weno3 and returns recompute vector                     --*/
/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
{

  Scalar_t *restrict LP, *restrict LN, *restrict alpha, *restrict flux;

  Scalar_t upwind = 1.0;

  Scalar_t weno_eps;
  Scalar_t D_C_CPt;
  Scalar_t D_C_MCt;
  Scalar_t D_M_Tt;
  Scalar_t D_CP_Tt;
  Scalar_t D_P_Tt;
  Scalar_t D_MC_Tt;
  Scalar_t vel_abs_max = 0.0;
  Scalar_t cce;

  Scalar_t p0m, p1m, p0p, p1p;
  Scalar_t B0m, B1m, B0p, B1p;
  Scalar_t w0m, w1m, w0p, w1p;
  Scalar_t wm_sum, wp_sum;
  Scalar_t OM0m, OM1m, OM0p, OM1p;
  Scalar_t um, up;
  Index_t species, energy, mu;

  D_C_CPt = 1/2.0;
  D_C_MCt = 1/2.0;
  D_M_Tt  = 1/3.0;
  D_CP_Tt = 2/3.0;
  D_P_Tt  = 1/3.0;
  D_MC_Tt = 2/3.0;
  weno_eps = 10.0*sqrt(DBL_MIN);

  flux  = (Scalar_t *) malloc(sizeof(Scalar_t) * NUM_SPECIES * NUM_ESTEPS * (NUM_MUSTEPS+1));
  LP    = (Scalar_t *) malloc(sizeof(Scalar_t) * NUM_SPECIES * NUM_ESTEPS * NUM_MUSTEPS);
  LN    = (Scalar_t *) malloc(sizeof(Scalar_t) * NUM_SPECIES * NUM_ESTEPS * NUM_MUSTEPS);
  alpha = (Scalar_t *) malloc(sizeof(Scalar_t) * NUM_SPECIES * NUM_ESTEPS * NUM_MUSTEPS);

/* ****** Get alpha ****** */

/* NOTE: v_avg is on extended "half" mesh.
   This means it has a length in energy of NUM_MUSTEPS+1
   Indicies:
   X    O    X    O    X    O    X    O    X    O    X
   0    0    1    1    2    2    3    3    NE-1 NE-1 NE
   (v_avg,flux on X), (alpha,LP/N,f,f1 on O)
*/

/* ***** Get inner points. */
  for (species = 0; species < NUM_SPECIES; species++) {
    for (energy = 0; energy < NUM_ESTEPS; energy++) {
      for (mu = 2; mu < NUM_MUSTEPS-2; mu++) {
        vel_abs_max =
            fabs(v_avg[idx_spemp1(species,energy,mu  )]);
        if (fabs(v_avg[idx_spemp1(species,energy,mu-1)]) > vel_abs_max) vel_abs_max = fabs(v_avg[idx_spemp1(species,energy,mu-1)]);
        if (fabs(v_avg[idx_spemp1(species,energy,mu-2)]) > vel_abs_max) vel_abs_max = fabs(v_avg[idx_spemp1(species,energy,mu-2)]);
        if (fabs(v_avg[idx_spemp1(species,energy,mu+1)]) > vel_abs_max) vel_abs_max = fabs(v_avg[idx_spemp1(species,energy,mu+1)]);
        if (fabs(v_avg[idx_spemp1(species,energy,mu+2)]) > vel_abs_max) vel_abs_max = fabs(v_avg[idx_spemp1(species,energy,mu+2)]);
        if (fabs(v_avg[idx_spemp1(species,energy,mu+3)]) > vel_abs_max) vel_abs_max = fabs(v_avg[idx_spemp1(species,energy,mu+3)]);

        alpha[idx_spem(species,energy,mu)] = vel_abs_max;
      }
    }
  }

  /* ***** Edge cases. */

  for (species = 0; species < NUM_SPECIES; species++) {
    for (energy = 0; energy < NUM_ESTEPS; energy++) {
      mu = 0;
      vel_abs_max =
          fabs(v_avg[idx_spemp1(species,energy,mu  )]);
      if (fabs(v_avg[idx_spemp1(species,energy,mu+1)]) > vel_abs_max) vel_abs_max = fabs(v_avg[idx_spemp1(species,energy,mu+1)]);
      if (fabs(v_avg[idx_spemp1(species,energy,mu+2)]) > vel_abs_max) vel_abs_max = fabs(v_avg[idx_spemp1(species,energy,mu+2)]);
      if (fabs(v_avg[idx_spemp1(species,energy,mu+3)]) > vel_abs_max) vel_abs_max = fabs(v_avg[idx_spemp1(species,energy,mu+3)]);

      alpha[idx_spem(species,energy,mu)] = vel_abs_max;

      mu = 1;
      vel_abs_max =
          fabs(v_avg[idx_spemp1(species,energy,mu  )]);
      if (fabs(v_avg[idx_spemp1(species,energy,mu-1)]) > vel_abs_max) vel_abs_max = fabs(v_avg[idx_spemp1(species,energy,mu-1)]);
      if (fabs(v_avg[idx_spemp1(species,energy,mu+1)]) > vel_abs_max) vel_abs_max = fabs(v_avg[idx_spemp1(species,energy,mu+1)]);
      if (fabs(v_avg[idx_spemp1(species,energy,mu+2)]) > vel_abs_max) vel_abs_max = fabs(v_avg[idx_spemp1(species,energy,mu+2)]);
      if (fabs(v_avg[idx_spemp1(species,energy,mu+3)]) > vel_abs_max) vel_abs_max = fabs(v_avg[idx_spemp1(species,energy,mu+3)]);

      alpha[idx_spem(species,energy,mu)] = vel_abs_max;

      mu = NUM_MUSTEPS-2;
      vel_abs_max =
          fabs(v_avg[idx_spemp1(species,energy,mu)]);
      if (fabs(v_avg[idx_spemp1(species,energy,mu+1)]) > vel_abs_max) vel_abs_max = fabs(v_avg[idx_spemp1(species,energy,mu+1)]);
      if (fabs(v_avg[idx_spemp1(species,energy,mu+2)]) > vel_abs_max) vel_abs_max = fabs(v_avg[idx_spemp1(species,energy,mu+2)]);
      if (fabs(v_avg[idx_spemp1(species,energy,mu-1)]) > vel_abs_max) vel_abs_max = fabs(v_avg[idx_spemp1(species,energy,mu-1)]);
      if (fabs(v_avg[idx_spemp1(species,energy,mu-2)]) > vel_abs_max) vel_abs_max = fabs(v_avg[idx_spemp1(species,energy,mu-2)]);

      alpha[idx_spem(species,energy,mu)] = vel_abs_max;

      mu = NUM_MUSTEPS-1;
      vel_abs_max =
          fabs(v_avg[idx_spemp1(species,energy,mu  )]);
      if (fabs(v_avg[idx_spemp1(species,energy,mu+1)]) > vel_abs_max) vel_abs_max = fabs(v_avg[idx_spemp1(species,energy,mu+1)]);
      if (fabs(v_avg[idx_spemp1(species,energy,mu-1)]) > vel_abs_max) vel_abs_max = fabs(v_avg[idx_spemp1(species,energy,mu-1)]);
      if (fabs(v_avg[idx_spemp1(species,energy,mu-2)]) > vel_abs_max) vel_abs_max = fabs(v_avg[idx_spemp1(species,energy,mu-2)]);

      alpha[idx_spem(species,energy,mu)] = vel_abs_max;
    }
  }

  for (species = 0; species < NUM_SPECIES; species++) {
    for (energy = 0; energy < NUM_ESTEPS; energy++) {
      for (mu = 0; mu < NUM_MUSTEPS; mu++) {

        LP[idx_spem(species,energy,mu)] = 0.5 * f[idx_spem(species,energy,mu)] *
          v_avg[idx_spemp1(species,energy,mu+1)] - alpha[idx_spem(species,energy,mu)];

        LN[idx_spem(species,energy,mu)] = 0.5 * f[idx_spem(species,energy,mu)] *
          v_avg[idx_spemp1(species,energy,mu)] + alpha[idx_spem(species,energy,mu)];

      }
    }
  }

  for (species = 0; species < NUM_SPECIES; species++) {
    for (energy = 0; energy < NUM_ESTEPS; energy++) {
      for (mu = 2; mu < NUM_MUSTEPS-1; mu++) {

        p0m = (1.0 + D_C_MCt)*LN[idx_spem(species,energy,mu-1)] - D_C_MCt*LN[idx_spem(species,energy,mu-2)];
        p1m =        D_C_MCt *LN[idx_spem(species,energy,mu-1)] + D_C_CPt*LN[idx_spem(species,energy,mu  )];
        p0p = (1.0 + D_C_CPt)*LP[idx_spem(species,energy,mu  )] - D_C_CPt*LP[idx_spem(species,energy,mu+1)];
        p1p =        D_C_CPt *LP[idx_spem(species,energy,mu  )] + D_C_MCt*LP[idx_spem(species,energy,mu-1)];

        B0m = 4.0*pow((D_C_MCt*(LN[idx_spem(species,energy,mu-1)] - LN[idx_spem(species,energy,mu-2)])),2);
        B1m = 4.0*pow((D_C_CPt*(LN[idx_spem(species,energy,mu  )] - LN[idx_spem(species,energy,mu-1)])),2);
        B0p = 4.0*pow((D_C_CPt*(LP[idx_spem(species,energy,mu+1)] - LP[idx_spem(species,energy,mu  )])),2);
        B1p = 4.0*pow((D_C_MCt*(LP[idx_spem(species,energy,mu  )] - LP[idx_spem(species,energy,mu-1)])),2);


        w0m = D_P_Tt /pow((weno_eps + B0m),2);
        w1m = D_MC_Tt/pow((weno_eps + B1m),2);
        w0p = D_M_Tt /pow((weno_eps + B0p),2);
        w1p = D_CP_Tt/pow((weno_eps + B1p),2);

        wm_sum = w0m + w1m;
        wp_sum = w0p + w1p;

        OM0m = w0m/wm_sum;
        OM1m = w1m/wm_sum;
        OM0p = w0p/wp_sum;
        OM1p = w1p/wp_sum;

        um = OM0m*p0m + OM1m*p1m;
        up = OM0p*p0p + OM1p*p1p;

        flux[idx_spemp1(species,energy,mu)] = up + um;

      }
    }
  }

  // Use upwinding for fluxes near the boundary.

  for (species = 0; species < NUM_SPECIES; species++) {
    for (energy = 0; energy < NUM_ESTEPS; energy++) {

      mu = 1;

      cce = copysign(upwind,v_avg[idx_spemp1(species,energy,mu)]);
      flux[idx_spemp1(species,energy,mu)] = v_avg[idx_spemp1(species,energy,mu)]
             * 0.5 * ( (1.0 - cce) * f[idx_spem(species,energy,mu  )]
                   +   (1.0 + cce) * f[idx_spem(species,energy,mu-1)]);

      mu = NUM_MUSTEPS-1;

      cce = copysign(upwind,v_avg[idx_spemp1(species,energy,mu)]);
      flux[idx_spemp1(species,energy,mu)] = v_avg[idx_spemp1(species,energy,mu)]
             * 0.5 * ( (1.0 - cce) * f[idx_spem(species,energy,mu  )]
                   +   (1.0 + cce) * f[idx_spem(species,energy,mu-1)]);

    }
  }

  for (species = 0; species < NUM_SPECIES; species++) {
    for (energy = 0; energy < NUM_ESTEPS; energy++) {
      for (mu = 1; mu < NUM_MUSTEPS-1; mu++) {

        f1[idx_spem(species,energy,mu)]=(flux[idx_spemp1(species,energy,mu+1)] -
                                         flux[idx_spemp1(species,energy,mu  )])/dmu;

      }
    }
  }

  // Boundary conditions:
  // Solid wall
  // For inflow, allow stuff to enter boundary cell, but not leave it.
  // For outflow, allow stuff to leave boundary, nothing enters.

  for (species = 0; species < NUM_SPECIES; species++) {
    for (energy = 0; energy < NUM_ESTEPS; energy++) {

      // Left boundary

      mu = 0;

      if (v_avg[idx_spemp1(species,energy,mu+1)] <= 0) {
        f1[idx_spem(species,energy,mu)] = v_avg[idx_spemp1(species,energy,mu+1)]
                                          *f[idx_spem(species,energy,mu+1)]/dmu;
      }
      else
      {
        f1[idx_spem(species,energy,mu)] = v_avg[idx_spemp1(species,energy,mu+1)]
                                          *f[idx_spem(species,energy,mu)]/dmu;
      }

      // Right boundary

      mu = NUM_MUSTEPS-1;

      if (v_avg[idx_spemp1(species,energy,mu)] >= 0) {
        f1[idx_spem(species,energy,mu)] = v_avg[idx_spemp1(species,energy,mu)]
                                          *(-f[idx_spem(species,energy,mu-1)])/dmu;
      }
      else
      {
        f1[idx_spem(species,energy,mu)] = v_avg[idx_spemp1(species,energy,mu)]
                                          *(-f[idx_spem(species,energy,mu)])/dmu;
      }
    }
  }

  free(flux);
  free(LP);
  free(LN);
  free(alpha);

} /*-------- END AdiabaticFocusing_Operator_WENO3( ) -----------------*/
/*------------------------------------------------------------------*/

/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/*--*/             void                                         /*--*/
/*--*/     GetStreamList( Index_t iterIndex )                   /*--*/
/*--                                                              --*/
/*-- streamlist with removed duplicates                           --*/
/*-- gathers stream at [face][row][col+mpi_rank] onto this proc   --*/
/*------------------------------------------------------------------*/
{/*-----------------------------------------------------------------*/

  Index_t shell, down, up, workIndex;
  Index_t shellP;
  Vec_t r0, r1;
  Scalar_t dsh;
  Scalar_t dshmin;

  workIndex = mpi_rank + N_PROCS * iterIndex;

  dshmin = config.dsh_min;

  if ((config.masHelCouple > 0) &&
      (t_global > masTime[config.masNumFiles-1])) {
      dshmin = config.dsh_hel_min;
  }

  // check if this process received work
  if (workIndex < NUM_STREAMS)
  {

    // initialize variables and arrays
    streamlistSize = 0;
    down = 0;

    for (shell = 0 ; shell < TOTAL_NUM_SHELLS;  shell++ ) {

      if (shell < (TOTAL_NUM_SHELLS - 1)) {

        up = shell + 1;

      } else {

        up = shell;
        down = shell - 1;

      }

      r1.x = streamGrid[up].r.x * config.rScale;
      r1.y = streamGrid[up].r.y * config.rScale;
      r1.z = streamGrid[up].r.z * config.rScale;

      r0.x = streamGrid[down].r.x * config.rScale;
      r0.y = streamGrid[down].r.y * config.rScale;
      r0.z = streamGrid[down].r.z * config.rScale;

      dsh = sqrt( (r1.x-r0.x)*(r1.x-r0.x) +
                  (r1.y-r0.y)*(r1.y-r0.y) +
                  (r1.z-r0.z)*(r1.z-r0.z)   );

      if (( dsh > dshmin ) || (shell == TOTAL_NUM_SHELLS - 1) ) {

        shellList[streamlistSize] = shell;
        streamlistSize++;
        shellRef[shell] = shell;
        down = shell + 1;

      } else {

        shellRef[shell] = shell + 1;

        shellP = shell - 1;

        while ( (shellP >= 0 ) && (shellRef[shellP] != shellP) ) {

          shellRef[shellP] = shellRef[shellP + 1];
          shellP--;

        }

      }

    }

//  Update ds variables with new node spacing.
    FindSegmentLengths();

  }

}/*-------- END GetStreamList() ------------------------------------*/
/*------------------------------------------------------------------*/

/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/*--*/             int                                          /*--*/
/*--*/     DiffuseStreamData( Index_t iterIndex,                /*--*/
/*--*/                        Scalar_t dt )                     /*--*/
/*--                                                              --*/
/*-- Calculate diffusion along a streamline                       --*/
/*------------------------------------------------------------------*/
{/*-----------------------------------------------------------------*/

  Index_t   shell, workIndex, face, row, col;
  Index_t   nsteps, step, species, energy, mu, slist;
  Scalar_t  NUM_MUSTEPS_I;
  Scalar_t  dtMin, dtProp,dtProp_i,vgrid_current,vgrid_current_i;
  Scalar_t  del_fac, iso, tau, rig;

  Scalar_t  *restrict f_old,          *restrict f_new;
  Scalar_t  *restrict iso_vec;
  Scalar_t  *restrict sep_seed_vec,   *restrict ds_i_multiplier_vec;
  Scalar_t  *restrict exp_mdt_tau_vec,*restrict del_fac_vec;

  const double one  = 1.0;

  NUM_MUSTEPS_I = one/NUM_MUSTEPS;

  workIndex = mpi_rank + N_PROCS * iterIndex;

  if (workIndex < NUM_STREAMS)
  {

    face = computeLines[workIndex][0];
    row  = computeLines[workIndex][1];
    col  = computeLines[workIndex][2];
//
// ****** Allocate temporary arrays ******
//
    f_old = malloc(streamlistSize*NUM_MUSTEPS*sizeof(Scalar_t));
    f_new = malloc(streamlistSize*NUM_MUSTEPS*sizeof(Scalar_t));

    sep_seed_vec         = malloc(3*sizeof(Scalar_t));
    ds_i_multiplier_vec  = malloc(streamlistSize*sizeof(Scalar_t));
    exp_mdt_tau_vec      = malloc(streamlistSize*sizeof(Scalar_t));
    iso_vec              = malloc(streamlistSize*sizeof(Scalar_t));
    del_fac_vec          = malloc(NUM_MUSTEPS*sizeof(Scalar_t));
//
// ****** Pre-load independent calculations invloving mu.
//

    for (species = 0; species < NUM_SPECIES; species++)
    {
      for (energy = 0; energy < NUM_ESTEPS; energy++)
      {
//
// ****** Get current particle velocity and rigidity^p.
//
        vgrid_current   =     vgrid[energy];
        vgrid_current_i = one/vgrid[energy];
        rig             =  rigidity[energy];
//
// ****** Compute stable time-step for explicit Euler sub-cycles.
//
        dtMin    = 0.4*dsMin*vgrid_current_i;


        // Need to modify V based on Nathans new document.
        // After modification, need to recalculate stable Euler dt...


        nsteps   = (Index_t) floor(dt/dtMin + one);
        dtProp   = dt/(one*nsteps);
        dtProp_i = one/dtProp;
//
// ****** Pre-load sub-cycle-step independent values.
//
        for (slist = 0; slist < streamlistSize;  slist++)
        {
          shell = shellList[slist];

// ****** Modifiy multiplier based on time-scale of mean-free-path.
          tau = vgrid_current_i*rig*pow(streamGrid[shell].rmag*config.rScale,
                                        config.mfpRadialPower)*config.lamo;
//          tau = vgrid_current_i*rig*config.lamo;
//          tau = vgrid_current_i*rig*(config.mhdBsAu/streamGrid[shell].mhdBmag)*config.lamo;
          if (dtProp > tau){
             ds_i_multiplier_vec[slist] = ds_i[slist]*tau*dtProp_i;
          }else{
             ds_i_multiplier_vec[slist] = ds_i[slist]*one;
          }

// ****** Pre-compute expensive operations used in isotropize.
          exp_mdt_tau_vec[slist] = exp(-dtProp/tau);
        }

// ****** Save initial seed population [for use with BCs].
        for ( slist = 0; slist < 3; slist++ ){
          shell = shellList[slist];
          sep_seed_vec[slist] = sepSeedFunction(egrid[idx_se(species,energy)],
                                                streamGrid[shell].rmag);
        }

        for (mu = 0; mu < NUM_MUSTEPS; mu++){
          del_fac_vec[mu] = vgrid_current*dtProp*mugrid[mu];
        }
//
// ****** Load stream data into stride-1 arrays (for speed).
//
        for (mu = 0; mu < NUM_MUSTEPS; mu++){
          for (slist = 0; slist < streamlistSize; slist++){
            f_old[streamlistSize*mu + slist] =
              ePartsStream[idx_sspem(shellList[slist],species,energy,mu)];
          }
        }
//
// *********************************************************************
// ****** Start integration of Euler sub-cycle.
// *********************************************************************
//
        for ( step = 0; step < nsteps; step++ ){
//
// ****** Set the boundary condition.
//
          if ( (config.useEPBoundary > 0) && (config.useBoundaryFunction > 0) ){
            for (mu = 0; mu < NUM_MUSTEPS; mu++){
              for (slist = 0 ; slist < 3;  slist++){
                f_old[streamlistSize*mu + slist] = sep_seed_vec[slist];
              }
            }
          }
//
// ****** Apply the upwinded advection step.
//
          for ( mu = 0; mu < NUM_MUSTEPS; mu++ ){
            del_fac = del_fac_vec[mu];

            if ( mugrid[mu] >= 0.0 ){
// ****** Unroll slist=0 loop iteration to allow vectorization.
              f_new[streamlistSize*mu] = f_old[streamlistSize*mu];

              for ( slist = 1; slist < streamlistSize; slist++ ){
                f_new[streamlistSize*mu + slist] =
                         f_old[streamlistSize*mu + slist]
                       + ds_i_multiplier_vec[slist] * del_fac
                       * ( f_old[streamlistSize*mu + slist - 1]
                         - f_old[streamlistSize*mu + slist] );
              }
            }
            else{
// ****** Unroll slist=(streamlistSize-1) loop iteration to allow vectorization.
              f_new[streamlistSize*mu + (streamlistSize-1)] =
              f_old[streamlistSize*mu + (streamlistSize-1)];

              for ( slist = 0; slist < (streamlistSize-1); slist++){
                f_new[streamlistSize*mu + slist] =
                         f_old[streamlistSize*mu + slist]
                       + ds_i_multiplier_vec[slist] * del_fac
                       * ( f_old[streamlistSize*mu + slist]
                         - f_old[streamlistSize*mu + slist + 1] );
              }
            }
          }
//
// ****** Isotropize across mu levels.
//
          for ( slist = 0 ; slist < streamlistSize;  slist++ ){
            iso_vec[slist] = 0.0;
          }

          for ( mu = 0; mu < NUM_MUSTEPS; mu++ ){
            for ( slist = 0 ; slist < streamlistSize;  slist++ ){
              iso_vec[slist] = iso_vec[slist] + f_new[streamlistSize*mu + slist];
            }
          }

          for ( mu = 0; mu < NUM_MUSTEPS; mu++ ){
            for ( slist = 0 ; slist < streamlistSize;  slist++ ){
              iso = NUM_MUSTEPS_I*iso_vec[slist];

// ****** Include resetting array in this last calculation.
              f_old[streamlistSize*mu + slist] = iso +
               (f_new[streamlistSize*mu + slist] - iso) * exp_mdt_tau_vec[slist];
            }
          }

          // Should the bondary conditions above be applied here too?


        } /*-- END OF SUB-CYCLE STEPS --*/
// *********************************************************************

//
// ****** Load stream data back into eprem structure.
// ****** This uses f_old since it has already been reset.
//
        for (mu = 0; mu < NUM_MUSTEPS; mu++)
        {
          for (slist = 0; slist < streamlistSize; slist++)
          {
            ePartsStream[idx_sspem(shellList[slist],species,energy,mu)]
              = f_old[streamlistSize*mu + slist];
          }
        }
//
// ****** Assign non-unique components and check for NaNs,Infs, and negative values.
//
        for ( mu = 0; mu < NUM_MUSTEPS; mu++ )
        {
          for ( shell = 0; shell < TOTAL_NUM_SHELLS; shell++ )
          {
            if (shell != shellRef[shell]){
              ePartsStream[idx_sspem(shell,species,energy,mu)] =
                ePartsStream[idx_sspem(shellRef[shell],species,energy,mu)];
            }
            checkNaN(mpi_rank, face, row, col, shell,
                     ePartsStream[idx_sspem(shell,species,energy,mu)],
                     "Diffuse Stream Data");
            checkInf(mpi_rank, face, row, col, shell,
                     ePartsStream[idx_sspem(shell,species,energy,mu)],
                     "Diffuse Stream Data");
            if ( ePartsStream[idx_sspem(shell,species,energy,mu)] < DBL_MIN )
            {
            //  warn(face, row, col, shell, species, energy, mu,
            //       "DiffuseStreamData: distribution less than DBL_MIN",
            //       &ePartsStream[idx_sspem(shell,species,energy,mu)]);
              ePartsStream[idx_sspem(shell,species,energy,mu)] = DBL_MIN;
            }

          }
        }

      } /*-- ENERGY --*/
    } /*-- SPECIES --*/
//
// ****** Free temporary array memory.
//
    free(f_old);
    free(f_new);
    free(sep_seed_vec);
    free(ds_i_multiplier_vec);
    free(exp_mdt_tau_vec);
    free(del_fac_vec);
    free(iso_vec);

  }

  return 0;

}/*-------- END DiffuseStreamData() --------------------------------*/


/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/*--*/     void                                                 /*--*/
/*--*/     FindSegmentLengths( void )                           /*--*/
/*--                                                              --*/
/*-- find array of displacements between nodes                    --*/
/*------------------------------------------------------------------*/
{/*-----------------------------------------------------------------*/

  Vec_t r0, r1, r9;
  Scalar_t dsh, ds1, ds9;
  Index_t thisShell, nextShell, lastShell;
  Index_t slist;
  const double one = 1.0;
  const double half = 0.5;

  dsMin = 10000.0 * config.rScale;

  for (slist = 0 ; slist < streamlistSize;  slist++ ) {

    thisShell = shellList[slist];
    /*     printf("thisShell = %d\n",thisShell); */

    if (slist == 0 ) {

      nextShell = shellList[slist + 1];

      /*  printf("nextShell = %d\n",nextShell); */

      r1.x = streamGrid[nextShell].r.x * config.rScale;
      r1.y = streamGrid[nextShell].r.y * config.rScale;
      r1.z = streamGrid[nextShell].r.z * config.rScale;

      r0.x = streamGrid[thisShell].r.x * config.rScale;
      r0.y = streamGrid[thisShell].r.y * config.rScale;
      r0.z = streamGrid[thisShell].r.z * config.rScale;

      dsh = sqrt((r1.x - r0.x) * (r1.x - r0.x) +
                 (r1.y - r0.y) * (r1.y - r0.y) +
                 (r1.z - r0.z) * (r1.z - r0.z));

      if (dsh < dsMin) {
        dsMin = dsh;
      }

      ds[slist] = dsh;

    } else if (slist == ( streamlistSize - 1 )) {

      lastShell = shellList[slist-1];
      /* printf("lastShell = %d\n",lastShell); */

      r9.x = streamGrid[lastShell].r.x * config.rScale;
      r9.y = streamGrid[lastShell].r.y * config.rScale;
      r9.z = streamGrid[lastShell].r.z * config.rScale;

      r0.x = streamGrid[thisShell].r.x * config.rScale;
      r0.y = streamGrid[thisShell].r.y * config.rScale;
      r0.z = streamGrid[thisShell].r.z * config.rScale;

      dsh = sqrt((r0.x - r9.x) * (r0.x - r9.x) +
                 (r0.y - r9.y) * (r0.y - r9.y) +
                 (r0.z - r9.z) * (r0.z - r9.z));

      if (dsh < dsMin) {
        dsMin = dsh;
      }

      ds[slist] = dsh;

    } else {

      lastShell = shellList[slist-1];
      nextShell = shellList[slist+1];

      r1.x = streamGrid[nextShell].r.x * config.rScale;
      r1.y = streamGrid[nextShell].r.y * config.rScale;
      r1.z = streamGrid[nextShell].r.z * config.rScale;

      r0.x = streamGrid[thisShell].r.x * config.rScale;
      r0.y = streamGrid[thisShell].r.y * config.rScale;
      r0.z = streamGrid[thisShell].r.z * config.rScale;

      r9.x = streamGrid[lastShell].r.x * config.rScale;
      r9.y = streamGrid[lastShell].r.y * config.rScale;
      r9.z = streamGrid[lastShell].r.z * config.rScale;

      ds1 = sqrt((r1.x - r0.x) * (r1.x - r0.x) +
                 (r1.y - r0.y) * (r1.y - r0.y) +
                 (r1.z - r0.z) * (r1.z - r0.z));

      ds9 = sqrt((r0.x - r9.x) * (r0.x - r9.x) +
                 (r0.y - r9.y) * (r0.y - r9.y) +
                 (r0.z - r9.z) * (r0.z - r9.z));

      if (ds1 < dsMin) {
        dsMin = ds1;
      }
      if (ds9 < dsMin) {
        dsMin = ds9;
      }
      ds[slist] = half * (ds1 + ds9);

    } /*-- endif (shell) --*/

    if (dsMin <= 0.0)
      panic("FindSegmentLengths(): Underflow in ds");

    ds_i[slist] = one / ds[slist];

  } /* shell loop */

}/*-------- END FindSementLengths() --------------------------------*/
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/*--*/             void                                         /*--*/
/*--*/     isotropize(Scalar_t dt,
                      Index_t species,
                      Index_t energy,
                      Index_t shell )
/*--                                                              --*/
/*--  returns a relaxed distribution over the relaxation time     --*/
/*------------------------------------------------------------------*/
{/*-----------------------------------------------------------------*/

  Scalar_t tau;
  Scalar_t iso;
  Index_t mu;

  // r = streamGrid[shell].rmag * config.rScale;

  iso = 0.0;
  tau = 0.0;

  for (mu = 0; mu < NUM_MUSTEPS; mu++)
    iso += ePartsStream[idx_sspem(shell, species, energy, mu)];

  iso /= (1.0 * NUM_MUSTEPS);

  tau = meanFreePath(species, energy,
                     streamGrid[shell].rmag*config.rScale)/vgrid[energy];

  for (mu = 0; mu < NUM_MUSTEPS; mu++)
    ePartsStream[idx_sspem(shell, species, energy, mu)] =
        iso + (ePartsStream[idx_sspem(shell, species, energy, mu)] - iso) * exp(-dt / tau);

} /*-------- END isotropize( ) -------------------------------------*/
/*------------------------------------------------------------------*/

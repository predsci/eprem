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

#include <stdio.h>
#include <math.h>
#include <float.h>
#include "global.h"
#include "configuration.h"
#include "energeticParticles.h"
#include "energeticParticlesBoundary.h"
#include "cubeShellStruct.h"
#include "simCore.h"
#include "flow.h"
#include "error.h"

/*-----------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------*/
/*--*/     void                                                            /*--*/
/*--*/     seedPopulationCheck( Index_t iterIndex )                        /*--*/
/*-----------------------------------------------------------------------------*/
{/*----------------------------------------------------------------------------*/

  Index_t workIndex, shell, species, energy, mu;

  Scalar_t distFunction;

  workIndex = mpi_rank + N_PROCS * iterIndex;

  if (workIndex < NUM_STREAMS)
  {

    for (shell = 0; shell < TOTAL_NUM_SHELLS; shell++) {
      for (species = 0; species < NUM_SPECIES; species++) {
        for (energy = 0; energy < NUM_ESTEPS; energy++) {

          distFunction = sepSeedFunction( egrid[energy], streamGrid[shell].rmag);

          for (mu = 0; mu < NUM_MUSTEPS; mu++)
            if (ePartsStream[idx_sspem(shell,species,energy,mu)] < distFunction)
              ePartsStream[idx_sspem(shell,species,energy,mu)] = distFunction;

        }

      }

    }

  }

}
/*--------------- end seedPopulationCheck() -----------------------------------*/
/*-----------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------*/
/*--*/     Scalar_t                                                        /*--*/
/*--*/     sepSeedFunction( Scalar_t energy, Scalar_t r)                   /*--*/
/*-----------------------------------------------------------------------------*/
{/*----------------------------------------------------------------------------*/

  Scalar_t normJHe, normJ;
  Scalar_t normJ0;
  Scalar_t normRadius, radialTerm;
  Scalar_t normE1, powerLawTerm;
  Scalar_t normE0, expTerm;
  Scalar_t normf;

  const double two = 2.0;
  const double r1  = 1.0; // Reference radius = 1 au
  const double E1  = 1.0; // Reference energy = 1 MeV/nuc
  const double efac = MEV / (MP * C * C);

  Scalar_t J0    = config.boundaryFunctAmplitude;
  Scalar_t xi    = config.boundaryFunctXi;
  Scalar_t gamma = config.boundaryFunctGamma;
  Scalar_t beta  = config.boundaryFunctBeta;
  Scalar_t E0    = config.boundaryFunctEcutoff;

  //if (energy/efac < 1.0){
  //  energy = 1.0 * efac;
  //}

  normJ0 = J0 * (MP * C) / (MHD_DENSITY_NORM * MEV);
  normRadius = r1 / config.rScale;
  normE1 = E1 * efac;
  normE0 = E0 * efac;

  radialTerm   = pow( (r / normRadius), -beta );
  powerLawTerm = pow( (energy / normE1), -gamma);
  expTerm      = exp( -(energy / normE0) );

  normJHe = normJ0 * radialTerm * powerLawTerm * expTerm;
  normJ  = normJHe / xi;

  normf = normJ / (two * energy);

  // If value is lower than the minimum double precision value
  // set it to that value (avoid 0).
  // This allows the log to be taken in adiabatic change.

  if (normf < DBL_MIN) normf = DBL_MIN;

  return normf;

}
/*-----------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------*/

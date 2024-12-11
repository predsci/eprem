/*-----------------------------------------------
-- EMMREM: cubeShellInit.h
--
-- Grid data structure initialization.
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

#ifndef CUBESHELLINIT_H
#define CUBESHELLINIT_H

#include "baseTypes.h"

void gridStructInit( void ); /*-- Init. the grid.             --*/
void initNEWS( Index_t shell ); /*-- Grid cells get NEWS links.  --*/
void checkNEWS( Index_t shell ); /*-- Grid cells check NEWS links.  --*/
void initCubeCoords( Index_t shell ); /*-- Set unit-cube coords.       --*/
void initSphereCoordsFromOuterBoundary( Index_t shell ); /*-- Set unit-cube coords.       --*/
void initSphereCoords_user( Index_t shell ); /*-- Scale to unit sphere.       --*/
void initSphereCoords( Index_t shell ); /*-- Scale to unit sphere.       --*/
void initCopyAll( void ); /*-- Copy inner shell to all.    --*/
void initStreamNeighbors(void ); /*-- Set all stream links.       --*/
void initInnerShellValues( void );      /*-- Initializing the values on the inner-most shell --*/

#endif

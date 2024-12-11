/*-----------------------------------------------
-- EMMREM: baseTypes.h
--
-- General base types and other general usage items.
-- MPI types are defined based on the types appearing here. Any modifications
-- to these types require corresponding modifications to MPI type
-- initialization routines in baseTypes.c.
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

#ifndef BASETYPES_H
#define BASETYPES_H

#include <mpi.h>

typedef double            doublereal;   /*-- See f2c.h.--*/
typedef long int          integer;      /*-- See f2c.h.--*/
//typedef long int          Index_t;
typedef int               Index_t;
typedef double            Time_t;
typedef double            Coord_t;
typedef double            Scalar_t;
typedef double            Radian_t;
typedef int               Bool_t;
typedef int               Flag_t;
typedef int               MPI_Rank_t;
typedef int               MPI_Flag_t;

#define T 1
#define F 0

/*-- [see numpy.pi] --*/
#define PI  3.1415926535897932
#define TWO_PI 6.283185307179586

#define VERYSMALL 1.0e-33
#define BADVALUE 1.0e33
#define BADINT 2147483647

#define MP  1.6726e-24
      /*-- grams                                 --*/
#define EV  1.6022e-12
      /*-- ergs                                  --*/
#define MEV 1.6022e-6
      /*-- ergs                                  --*/
#define GEV 1.6022e-3
#define Q   4.80320425e-10
      /*-- electron charge [statcoul]            --*/
#define C   2.99792458e10
      /*-- speed of light  [cm/s    ]            --*/
#define MZERO (MP*C*C/GEV)
/*-- mp*c^2 in GeV                         --*/
#define AU 1.495978707e13
      /*-- 1 Astronomical Unit [cm  ]            --*/

#define RSUN 6.96e10
     /*-- Solar radius [cm  ]            --*/

#define RSAU (RSUN/AU)
    /*-- Solar radius [in AU]            
      -- 0.00464913034353770 [AU] */
    
#define TAU ( AU / C )
      /*-- Simulation Timescale in s             --
        -- 1 timestep = 1 TAU                    --
        -- TAU = AU/c =   8.31674639726927 minutes --
        --            = 499.00478383615643 seconds --
*/
        
#define DAY ( TAU / (24.0 * 60.0 * 60.0)  )
      /* -- 1 timestep*DAY = Julian Days --
         -- DAY = 0.005775518331436995 -- */

#define MHD_DENSITY_NORM 1.0
      /*-- Density norm in cm-3                  --*/
#define MHD_B_NORM ( sqrt(MP*MHD_DENSITY_NORM) * C )
      /*-- Field strength at 1 AU                --*/ 
      /*-- Bnorm^2 = rho0 c^2                    --*/
      /*-- (vA/c)^2= (B/Bnorm)^2/(4 pi)          --*/
      /*-- rho0 = 1.0 /cm^3                      --*/
      /*-- 1.29228e-12 = sqrt(mp)                --*/

/* normalization for ion gyrofrequency */
#define OM ( MHD_B_NORM * Q * AU / (MP * C * C) )

#define FCONVERT ( 1.0e30 / (C * C * C * C) )
/* used to convert to f from j into s^3/km^6 --*/

#define  VOLT 0.33333e-2
/* 1 VOLT = 0.333e-2 statvolt , pc/q in statvolt */

#define THRESH 0.025
/* threshold for perp diff & drift shell */

#define MAS_TIME_NORM 1445.87003080685
#define MAS_LENGTH_NORM 6.96e10
#define MAS_RHO_NORM 1.6726e-16
#define MAS_TIME_CONVERT MAS_TIME_NORM / (24.0 * 60.0 * 60.0)
#define MAS_V_CONVERT ( MAS_LENGTH_NORM / MAS_TIME_NORM ) / C
#define MAS_RHO_CONVERT (MAS_RHO_NORM / MP) / MHD_DENSITY_NORM
#define MAS_B_CONVERT sqrt(4.0*PI*MAS_RHO_NORM)*MAS_LENGTH_NORM / MAS_TIME_NORM / MHD_B_NORM

#define MAX_STRING_SIZE 240

#define MHD_DEFAULT 0
#define MHD_ENLIL   1
#define MHD_LFMH    2
#define MHD_BATSRUS 3
#define MHD_MAS     4

/*.............................Vec_t...*/
typedef struct {
    Coord_t x;
    Coord_t y;
    Coord_t z;
}
/*.....................END..*/ Vec_t;

/*A spherical vector */
/*.....................SphVec_t...*/
typedef struct {
  Coord_t r;
  Coord_t theta;
  Coord_t phi;
}
/*.....................END..*/ SphVec_t; 

/*---*/         void                                   /*---*/
/*---*/   initMPI_baseTypes(void );                    /*---*/

extern  MPI_Datatype Index_T;
extern  MPI_Datatype Time_T;
extern  MPI_Datatype Coord_T;
extern  MPI_Datatype Scalar_T;
extern  MPI_Datatype Radian_T;
extern  MPI_Datatype Bool_T;
extern  MPI_Datatype Flag_T;
extern  MPI_Datatype Vec_T;
extern  MPI_Datatype SphVec_T;
extern  MPI_Datatype MPI_Rank_T;
extern  MPI_Datatype MPI_Flag_T;

#endif

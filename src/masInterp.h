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

#ifndef MASINTERP_H
#define MASINTERP_H

#define Bpindex(r,t,p) ((r) + masBprDimMax[0] * (t) + masBprDimMax[0] * masBptDimMax[0] * (p))
#define Btindex(r,t,p) ((r) + masBtrDimMax[0] * (t) + masBtrDimMax[0] * masBttDimMax[0] * (p))
#define Brindex(r,t,p) ((r) + masBrrDimMax[0] * (t) + masBrrDimMax[0] * masBrtDimMax[0] * (p))
#define Vpindex(r,t,p) ((r) + masVprDimMax[0] * (t) + masVprDimMax[0] * masVptDimMax[0] * (p))
#define Vtindex(r,t,p) ((r) + masVtrDimMax[0] * (t) + masVtrDimMax[0] * masVttDimMax[0] * (p))
#define Vrindex(r,t,p) ((r) + masVrrDimMax[0] * (t) + masVrrDimMax[0] * masVrtDimMax[0] * (p))
#define Dindex(r,t,p)  ((r) + masDrDimMax[0]  * (t) + masDrDimMax[0]  * masDtDimMax[0] * (p))

typedef struct {
  SphVec_t r;           /* Position in spherical coords. */
  SphVec_t mhdB;        /* B-field    */
  SphVec_t mhdV;        /* Velocity   */
  Scalar_t mhdD;        /* Density    */
  SphVec_t curlBoverB2; /* del x B/B^2 */
} masNode_t;

extern masNode_t masNode;
extern Index_t unwindPhiOffset;

int masBinarySearch(float *, float, int, int);

Scalar_t masTriLinearBinarySearch(float *, Scalar_t, int *, int *, int, int);

Scalar_t masInterpolate(float *,
                        int, int, int, int, int, int,
                        Scalar_t, Scalar_t, Scalar_t,
                        int, int);
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*--*/  SphVec_t                                                          /*--*/
/*--*/  fetchMasB( SphVec_t r,                                            /*--*/
/*--*/             float masBp0[], float masBt0[], float masBr0[],        /*--*/
/*--*/             float masBp1[], float masBt1[], float masBr1[],        /*--*/
/*--*/             Scalar_t s );                                          /*--*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*--*/ SphVec_t masCurlBoverB2( SphVec_t r,                               /*--*/
/*--*/                    float Bp0[], float Bt0[], float Br0[],          /*--*/
/*--*/                    float Bp1[], float Bt1[], float Br1[],          /*--*/
/*--*/                    int bp_r0, int bp_r1,                           /*--*/
/*--*/                    int bp_t0, int bp_t1,                           /*--*/
/*--*/                    int bp_p0, int bp_p1,                           /*--*/
/*--*/                    int bt_r0, int bt_r1,                           /*--*/
/*--*/                    int bt_t0, int bt_t1,                           /*--*/
/*--*/                    int bt_p0, int bt_p1,                           /*--*/
/*--*/                    int br_r0, int br_r1,                           /*--*/
/*--*/                    int br_t0, int br_t1,                           /*--*/
/*--*/                    int br_p0, int br_p1,                           /*--*/
/*--*/                    Scalar_t s);                                    /*--*/
/*--*/                                                                    /*--*/
/*--*/                                                                    /*--*/
/*--     Calculate the curl of B/B^2 for use in the drift velocity          --*/
/*--                                                                        --*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*--*/ void masTriLinear( SphVec_t position,                              /*--*/
/*--*/                    float *, float *, float *,                      /*--*/
/*--*/                    float *, float *, float *,                      /*--*/
/*--*/                    float *,                                        /*--*/
/*--*/                    float *, float *, float *,                      /*--*/
/*--*/                    float *, float *, float *,                      /*--*/
/*--*/                    float *,                                        /*--*/
/*--*/                    Scalar_t);                                      /*--*/
/*--*/                                                                    /*--*/
/*--*/                                                                    /*--*/
/*--  Get the data from the nearest mas node and interpolate                --*/
/*--  in time.                                                              --*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*--*/ void masHelTriLinear( SphVec_t position,                           /*--*/
/*--*/                    float *, float *, float *,                      /*--*/
/*--*/                    float *, float *, float *,                      /*--*/
/*--*/                    float *,                                        /*--*/
/*--*/                    float *, float *, float *,                      /*--*/
/*--*/                    float *, float *, float *,                      /*--*/
/*--*/                    float *,                                        /*--*/
/*--*/                    Scalar_t);                                      /*--*/
/*--*/                                                                    /*--*/
/*--*/                                                                    /*--*/
/*--  Get the data from the nearest mas node and interpolate                --*/
/*--  in time.                                                              --*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/

/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/*--*/ SphVec_t fieldAlignedFlow(SphVec_t B);                   /*--*/
/*--*/                                                          /*--*/
/*--   Field aligned flow                                         --*/
/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/*--*/ void																											/*--*/
/*--*/ masGetNode(SphVec_t position,                            /*--*/
/*--*/            Node_t node);                         /*--*/
/*--*/                                                          /*--*/
/*--   gets the mas data at the specified position                --*/
/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/

/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/*--*/ void masWind(Node_t node);                               /*--*/
/*--*/                                                          /*--*/
/*--  Revert to the Parker wind model.                            --*/
/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/

/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/*--*/ void rotateCoupledDomain( Scalar_t dt );                 /*--*/
/*--*/                                                          /*--*/
/*--  If inside the coupled domain, rotate the nodes              --*/
/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/

/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/*--*/ void resetDomainOffset( void );                          /*--*/
/*--*/                                                          /*--*/
/*--  If inside the coupled domain, rotate the nodes              --*/
/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/

/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
/*--*/    Vec_t                                                     /*---*/
/*--*/    vMas(Vec_t r, Node_t node );                              /*---*/
/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/

/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
/*--*/    Vec_t                                                     /*---*/
/*--*/    rungeKuttaFlow( Vec_t r0, Scalar_t dt, Node_t node );     /*---*/
/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/

/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
/*--*/    void                                                      /*---*/
/*--*/    masMoveNodes( Scalar_t dt );                              /*---*/
/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/


#endif

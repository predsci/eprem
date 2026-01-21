#ifndef CONFIGURATION_H
#define CONFIGURATION_H

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

#include <libconfig.h>
#include "cubeShellStruct.h"

typedef struct {

  Index_t  numNodesPerStream;
  Index_t  numRowsPerFace;
  Index_t  numColumnsPerFace;
  Index_t  numEnergySteps;
  Index_t  numMuSteps;
  Index_t  adiabaticChangeAlg;
  Index_t  adiabaticFocusAlg;

  Scalar_t  rScale;
  Scalar_t  flowMag;
  Scalar_t  mhdDensityAu;
  Scalar_t  mhdBAu;
  Scalar_t  simStartTime;
  Scalar_t  simStopTime;
  Scalar_t  tDel;
  Index_t   numEpSteps;
  Scalar_t  aziSunStart;
  Scalar_t  omegaSun;
  Scalar_t  lamo;
  Scalar_t  dsh_min;
  Scalar_t  dsh_hel_min;
  Index_t   mfpType;
  Scalar_t  mfpRadialPower;
  Scalar_t  mfpBPower;
  Scalar_t  rigidityPower;

  Scalar_t  kperxkpar;

  Scalar_t  eMin;
  Scalar_t  eMax;
  Index_t   useStochastic;
  Scalar_t  focusingLimit;
  Index_t   useEPBoundary;
  Index_t   checkSeedPopulation;

  Scalar_t  gammaEhigh;
  Scalar_t  gammaElow;

  Index_t   FailModeDump;

  Index_t   seedFunctionTest;

  Index_t   outputFloat;

  Index_t   unifiedOutput;
  Scalar_t  unifiedOutputTime;

  Index_t   pointObserverOutput;
  Scalar_t  pointObserverOutputTime;

  Index_t   streamFluxOutput;
  Scalar_t  streamFluxOutputTime;

  Index_t   subTimeCouple;

  Index_t    epremDomain;
  Scalar_t   epremDomainOutputTime;

  Index_t    unstructuredDomain;
  Scalar_t   unstructuredDomainOutputTime;

  Index_t    useAdiabaticChange;
  Scalar_t   dampCooling;
  Index_t    useAdiabaticFocus;
  Index_t    useShellDiffusion;
  Index_t    useParallelDiffusion;
  Index_t    useDrift;

  Index_t fluxLimiter;

  int useManualStreamSpawnLoc;
  Scalar_t* streamSpawnLocAzi;
  Scalar_t* streamSpawnLocZen;

  int numSpecies;
  Scalar_t*  mass;
  Scalar_t*  charge;

  int numObservers;
  Scalar_t*  obsR;
  Scalar_t*  obsTheta;
  Scalar_t*  obsPhi;

  Scalar_t idw_p;

  Index_t masTriLinear;
  Index_t masCouple;
  Index_t masCorRotateFake;
  Index_t masHelCouple;
  Index_t masNumFiles;
  Index_t masHelNumFiles;
  Index_t useMasSteadyStateDt;
  Index_t masSteadyState;
  char* masDirectory;
  char* masHelDirectory;
  Index_t masDigits;
  Index_t masHelDigits;

  Index_t masCoupledTime;
  Scalar_t masStartTime;
  Scalar_t epEquilibriumCalcDuration;

  Scalar_t preEruptionDuration;
  Index_t masInitMonteCarlo;
  Index_t masInitFromOuterBoundary;
  Scalar_t masInitRadius;
  Scalar_t masInitTimeStep;

  Scalar_t masRadialMin;
  Scalar_t masRadialMax;
  Scalar_t masHelRadialMin;
  Scalar_t masHelRadialMax;
  Scalar_t masVmin;

  Scalar_t parallelFlow;
  Index_t  fieldAligned;

  Scalar_t epCalcStartTime;

  Index_t masRotateSolution;

  Index_t useBoundaryFunction;
  Index_t boundaryFunctionInitDomain;

  Scalar_t boundaryFunctAmplitude;
  Scalar_t boundaryFunctXi;
  Scalar_t boundaryFunctGamma;
  Scalar_t boundaryFunctBeta;
  Scalar_t boundaryFunctEcutoff;

  Index_t   shockSolver;
  Scalar_t  shockDetectPercent;
  Scalar_t  minInjectionEnergy;
  Scalar_t  shockInjectionFactor;

  Index_t   idealShock;
  Scalar_t  idealShockSharpness;
  Scalar_t  idealShockScaleLength;
  Scalar_t  idealShockJump;
  Scalar_t  idealShockFalloff;
  Scalar_t  idealShockSpeed;
  Scalar_t  idealShockInitTime;
  Scalar_t  idealShockTheta;
  Scalar_t  idealShockPhi;
  Scalar_t  idealShockWidth;

  int dumpFreq;
  int outputRestart;
  int dumpOnAbort;
  int saveRestartFile;

  char * warningsFile;

  // these params are initialized from the above inputs
  Scalar_t   mhdUs;
  Scalar_t   mhdNsAu;
  Scalar_t   mhdBsAu;

  Scalar_t   simStartTimeDay;
  Scalar_t   simStopTimeDay;

} Config_t;

extern Config_t config;
extern config_t cfg;

void initGlobalParameters( char* configFilename );
void getParams( char* configFilename);
void setRuntimeConstants( void );

Index_t readInt(char *key, Index_t defaultVal, Index_t minVal, Index_t maxVal);
Scalar_t readDouble(char *key, Scalar_t defaultVal, Scalar_t minVal, Scalar_t maxVal);
const char *readString(char *key, char *defaultVal);
Scalar_t *readDoubleArray(char *key, int size, Scalar_t *defaultVal);

#endif

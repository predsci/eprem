/*-----------------------------------------------
 -- EMMREM: configuration.h
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <libconfig.h>
#include "global.h"
#include "configuration.h"
#include "mpiInit.h"
#include "error.h"

Config_t config;
config_t cfg;

const double third = 1.0/3.0;

void
initGlobalParameters( char* configFilename )
{

  getParams(configFilename);
  setRuntimeConstants();

}


void
getParams( char* configFilename)
{

  // initialize the config parser structure and read external config file
  config_init(&cfg);

  if(! config_read_file(&cfg, configFilename)) {

    if (mpi_rank == 0)
      printf("\n\nError on line %i of %s: %s\n\n", config_error_line(&cfg), configFilename, config_error_text(&cfg));

    panic("Unable to read configuration file.");

  }

  config.numNodesPerStream = readInt("numNodesPerStream",N_PROCS,N_PROCS,BADINT);
  config.numRowsPerFace = readInt("numRowsPerFace", 2, 1, BADINT);
  config.numColumnsPerFace = readInt("numColumnsPerFace", 2, 1, BADINT);
  config.numEnergySteps = readInt("numEnergySteps", 20, 2, BADINT);
  config.numMuSteps = readInt("numMuSteps", 20, 2, BADINT);

  config.rScale = readDouble("rScale", 0.005, VERYSMALL, BADVALUE);
  config.flowMag = readDouble("flowMag", 400.0e5, VERYSMALL, BADVALUE);
  config.mhdDensityAu = readDouble("mhdDensityAu", 8.30, VERYSMALL, BADVALUE);
  config.mhdBAu = readDouble("mhdBAu", 1.60e-5, VERYSMALL, BADVALUE);
  config.simStartTime = readDouble("simStartTime", 0.0, 0.0, BADVALUE);
  config.tDel = readDouble("tDel", 0.01041666666667, VERYSMALL, BADVALUE);
  config.simStopTime = readDouble("simStopTime", config.simStartTime + config.tDel, config.simStartTime, BADVALUE);
  config.numEpSteps = readInt("numEpSteps", 30, 1, BADINT);
  config.aziSunStart = readDouble("aziSunStart", 0.0, 0.0, BADVALUE);
  // This omega makes EPREM and MAS omega values equal to double precision.
  config.omegaSun = readDouble("omegaSun", (0.004144/MAS_TIME_NORM)*TAU, 0.0, BADVALUE);
  config.dsh_min = readDouble("dsh_min", 5.0e-5, VERYSMALL, BADVALUE);
  config.dsh_hel_min = readDouble("dsh_hel_min", 2.5e-4, VERYSMALL, BADVALUE);
  config.kperxkpar = readDouble("kperxkpar", 0.01, VERYSMALL, BADVALUE);
  config.lamo = readDouble("lamo", 1.0, VERYSMALL, BADVALUE);
  // Mean free path options:
  config.mfpType        = readInt("mfpType", 1, 1, 2);
  // ***  mfpType:  1 = radial scaling  2 = B-dependent.
  config.mfpRadialPower = readDouble("mfpRadialPower",2.0,VERYSMALL, BADVALUE);
  config.mfpBPower      = readDouble("mfpBPower",0.75,VERYSMALL, BADVALUE);
  config.rigidityPower = readDouble("rigidityPower", third, VERYSMALL, BADVALUE);
  config.focusingLimit = readDouble("focusingLimit", 1.0, 0.0, 1.0);

  config.eMin = readDouble("eMin", 1.0, VERYSMALL, BADVALUE);
  config.eMax = readDouble("eMax", 1000.0, config.eMin, BADVALUE);
  config.useStochastic = readInt("useStochastic", 0, 0, 1);
  config.useEPBoundary = readInt("useEPBoundary", 1, 0, 1);
  config.checkSeedPopulation = readInt("checkSeedPopulation", 1, 0, 1);

  config.seedFunctionTest = readInt("seedFunctionTest", 0, 0, 1);

  config.fluxLimiter = readInt("fluxLimiter", 1, 0, 1);

  config.gammaEhigh = readDouble("gammaEhigh", 0.0, -1.0 * BADVALUE, BADVALUE);
  config.gammaElow = readDouble("gammaElow", 0.0, -1.0 * BADVALUE, BADVALUE);

  config.FailModeDump = readInt("FailModeDump", 1, 0, 1);

  config.outputFloat = readInt("outputFloat", 0, 0, 1);

  config.unifiedOutput = readInt("unifiedOutput", 1, 0, 1);
  config.unifiedOutputTime = readDouble("unifiedOutputTime", 0.0, 0.0, BADVALUE);

  config.pointObserverOutput = readInt("pointObserverOutput", 0, 0, 1);
  config.pointObserverOutputTime = readDouble("pointObserverOutputTime", 0.0, 0.0, BADVALUE);

  config.streamFluxOutput = readInt("streamFluxOutput", 0, 0, 1);
  config.streamFluxOutputTime = readDouble("streamFluxOutputTime", 0.0, 0.0, BADVALUE);

  config.subTimeCouple = readInt("subTimeCouple", 0, 0, 1);

  config.epremDomain = readInt("epremDomain", 0, 0, 1);
  config.epremDomainOutputTime = readDouble("epremDomainOutputTime", 0.0, 0.0, BADVALUE);

  config.unstructuredDomain = readInt("unstructuredDomain", 0, 0, 1);
  config.unstructuredDomainOutputTime = readDouble("unstructuredDomainOutputTime", 0.0, 0.0, BADVALUE);

  config.useAdiabaticChange = readInt("useAdiabaticChange", 1, 0, 1);
  config.dampCooling = readDouble("dampCooling", 0.0, 0.0, 1.0);
  config.useAdiabaticFocus = readInt("useAdiabaticFocus", 1, 0, 1);
  config.useShellDiffusion = readInt("useShellDiffusion", 0, 0, 1);
  config.useParallelDiffusion = readInt("useParallelDiffusion", 1, 0, 1);
  config.useDrift = readInt("useDrift", 0, 0, 1);

  config.numSpecies = readInt("numSpecies", 1, 0, 100);
  Scalar_t defaultMass[1] = {1.0};
  Scalar_t defaultEnergy[1] = {1.0};
  config.mass = readDoubleArray("mass", config.numSpecies, defaultMass);
  config.charge = readDoubleArray("charge", config.numSpecies, defaultEnergy);

  config.numObservers = readInt("numObservers", 0, 0, 1000);
  Scalar_t defaultObserver[1] = {0};
  config.obsR = readDoubleArray("obsR", config.numObservers, defaultObserver);
  config.obsTheta = readDoubleArray("obsTheta", config.numObservers, defaultObserver);
  config.obsPhi = readDoubleArray("obsPhi", config.numObservers, defaultObserver);

  config.idw_p = readDouble("idw_p", 3.0, VERYSMALL, BADVALUE);

  config.masTriLinear = readInt("masTriLinear", 1, 0, 1);
  config.masCouple = readInt("masCouple", 0, 0, 1);
  config.masCorRotateFake = readInt("masCorRotateFake", 0, 0, 1);
  config.masHelCouple = readInt("masHelCouple", 0, 0, 1);
  config.masNumFiles = readInt("masNumFiles", 0, 0, 32767);
  config.masHelNumFiles = readInt("masHelNumFiles", 0, 0, 32767);
  config.useMasSteadyStateDt = readInt("useMasSteadyStateDt", 1, 0, 1);
  config.masSteadyState = readInt("masSteadyState", 1, 0, 1);
  config.masDirectory = (char*)readString("masDirectory"," ");
  config.masHelDirectory = (char*)readString("masHelDirectory"," ");
  config.masDigits = readInt("masDigits", 3, 0, 32767);
  config.masHelDigits = readInt("masHelDigits", 3, 0, 32767);

  config.masCoupledTime = readInt("masCoupledTime", 1, 0, 1);
  config.masStartTime = readDouble("masStartTime", 0.0, 0.0, BADVALUE);
  config.epEquilibriumCalcDuration = readDouble("epEquilibriumCalcDuration", 0.0, 0.0, BADVALUE);
  config.preEruptionDuration = readDouble("preEruptionDuration", 0.0, 0.0, BADVALUE);
  config.masRadialMin = readDouble("masRadialMin", 0.0, 0.0, BADVALUE);
  config.masRadialMax = readDouble("masRadialMax", 0.0, 0.0, BADVALUE);
  config.masHelRadialMin = readDouble("masHelRadialMin", 0.0, 0.0, BADVALUE);
  config.masHelRadialMax = readDouble("masHelRadialMax", 0.0, 0.0, BADVALUE);
  config.masVmin = readDouble("masVmin", 50.0e5, 0.0, BADVALUE);

  config.masInitFromOuterBoundary = readInt("masInitFromOuterBoundary", 2, 0, 2);
  config.masInitMonteCarlo = readInt("masInitMonteCarlo", 0, 0, 1);
  config.masInitRadius = readDouble("masInitRadius", 0.0, 0.0, BADVALUE);
  config.masInitTimeStep = readDouble("masInitTimeStep", 0.000011574074074, 0.0, BADVALUE);
  config.useManualStreamSpawnLoc = readInt("useManualStreamSpawnLoc",0,0,1);
  Scalar_t defaultPos[1] = {0.0};
  if (config.useManualStreamSpawnLoc > 0){
    config.streamSpawnLocAzi = readDoubleArray("streamSpawnLocAzi", 6*config.numRowsPerFace*config.numColumnsPerFace, defaultPos);
    config.streamSpawnLocZen = readDoubleArray("streamSpawnLocZen", 6*config.numRowsPerFace*config.numColumnsPerFace, defaultPos);
  }

  config.parallelFlow = readDouble("parallelFlow", 0.0, 0.0, BADVALUE);
  config.fieldAligned = readInt("fieldAligned", 0, 0, 1);

  config.epCalcStartTime = readDouble("epCalcStartTime", config.simStartTime, 0.0, BADVALUE);

  config.masRotateSolution = readInt("masRotateSolution", 1, 0, 1);

  config.useBoundaryFunction = readInt("useBoundaryFunction", 1, 0, 1);
  config.boundaryFunctionInitDomain = readInt("boundaryFunctionInitDomain", 1, 0, 1);

  config.boundaryFunctAmplitude = readDouble("boundaryFunctAmplitude", 1.0, VERYSMALL, BADVALUE);
  config.boundaryFunctXi = readDouble("boundaryFunctXi", 1.0, 0.0, BADVALUE);
  config.boundaryFunctGamma = readDouble("boundaryFunctGamma", 2.0, 0.0, BADVALUE);
  config.boundaryFunctBeta = readDouble("boundaryFunctBeta", 1.7, 0.0, BADVALUE);
  config.boundaryFunctEcutoff = readDouble("boundaryFunctEcutoff", 1.0, 0.0, BADVALUE);
  config.boundaryFunctR0 = readDouble("boundaryFunctR0", 1.0, config.rScale, BADVALUE);

  config.shockSolver = readInt("shockSolver", 0, 0, 1);
  config.shockDetectPercent = readDouble("shockDetectPercent", 1.0, 0.0, BADVALUE);
  config.minInjectionEnergy = readDouble("minInjectionEnergy", 0.01, VERYSMALL, BADVALUE);
  config.shockInjectionFactor = readDouble("shockInjectionFactor", 1.0, 0.0, BADVALUE);

  config.idealShock = readInt("idealShock", 0, 0, 1);
  config.idealShockSharpness = readDouble("idealShockSharpness", 1.0, VERYSMALL, BADVALUE);
  config.idealShockScaleLength = readDouble("idealShockScaleLength", 0.0046491, VERYSMALL, BADVALUE);
  config.idealShockJump = readDouble("idealShockJump", 4.0, VERYSMALL, BADVALUE);
  config.idealShockFalloff = readDouble("idealShockFalloff", 0.0, 0.0, BADVALUE);
  config.idealShockSpeed = readDouble("idealShockSpeed", 1500e5, VERYSMALL, BADVALUE);
  config.idealShockInitTime = readDouble("idealShockInitTime", config.simStartTime, config.simStartTime, BADVALUE);
  config.idealShockTheta = readDouble("idealShockTheta", 1.570796, 0.0, PI);
  config.idealShockPhi = readDouble("idealShockPhi", 0.0, 0.0, 2.0 * PI);
  config.idealShockWidth = readDouble("idealShockWidth", 0.0, 0.0, PI);

  config.dumpFreq = readInt("dumpFreq",1, 0, 1000000);
  config.outputRestart = readInt("outputRestart",0, 0, 1000000);
  config.dumpOnAbort = readInt("dumpOnAbort",0, 0, 1);
  config.saveRestartFile = readInt("saveRestartFile",0, 0, 1);

  config.warningsFile = (char*)readString("warningsFile", "warningsXXX.txt");

  config.adiabaticChangeAlg = readInt("adiabaticChangeAlg", 1, 1, 3);
  config.adiabaticFocusAlg = readInt("adiabaticFocusAlg", 1, 1, 3);

}


Index_t readInt(char *key, Index_t defaultVal, Index_t minVal, Index_t maxVal) {

  Index_t val;

  if (! config_lookup_int(&cfg, key, &val) )
    val = defaultVal;

  if ( (val < minVal) || (val > maxVal) ) {

    printf("%s is out of the acceptable range: %.4i <= %.4i <= %.4i\n", key, minVal, (Index_t)val, maxVal);
    panic("the configuration reader detected an invalid value.\n");

  }

  if (mpi_rank == 0)
    printf("%s: %i\n", key, (Index_t)val);

  return (Index_t)val;

}


Scalar_t readDouble(char *key, Scalar_t defaultVal, Scalar_t minVal, Scalar_t maxVal) {

  Scalar_t val;

  if (! config_lookup_float(&cfg, key, &val) )
    val = defaultVal;

  if ( (val < minVal) || (val > maxVal) ) {

    printf("%s is out of the acceptable range: %.4e < %.4e < %.4e\n", key, minVal, (Scalar_t)val, maxVal);
    panic("the configuration reader detected an invalid value.\n");

  }

  if (mpi_rank == 0)
    printf("%s: %.4e\n", key, (Scalar_t)val);

  return (Scalar_t)val;

}


const char *readString(char *key, char *defaultVal) {

  const char *val;

  if (! config_lookup_string(&cfg, key, &val) )
    val = defaultVal;

  if (mpi_rank == 0)
    printf("%s: %s\n", key, val);

  return val;

}


Scalar_t *readDoubleArray(char *key, int size, Scalar_t *defaultVal) {

  Index_t i;
  Scalar_t *val;
  const config_setting_t *Arr;

  val = (Scalar_t *)malloc(sizeof(double) * size);

  if (mpi_rank == 0)
    printf("%s: ", key);

  if (size > 0) {

    Arr = config_lookup(&cfg, key);

    for (i = 0; i < size; i++) {
      val[i] = config_setting_get_float_elem(Arr, i);
      if (mpi_rank == 0) printf("%.4e ", val[i]);
    }

    if (mpi_rank == 0) printf("\n");

    return val;

  } else {

    if (mpi_rank == 0){
//   RMC: THIS IS WRONG AND DOES NOT FIND THE LENGTH OF THE ARRAY!!!
//   NEED TO PASS THE DEFAULT SIZE IN AS WELL!!!
//      for (i = 0; i < (Index_t)(sizeof(defaultVal) / sizeof(*defaultVal)); i++)
//        printf("%.4e ", defaultVal[i]);
      printf("Using default values.  First element: %.4e ", defaultVal[0]);
    }

    if (mpi_rank == 0) printf("\n");

    return defaultVal;

  }

}


void
setRuntimeConstants( void )
{

  FACE_ROWS = config.numRowsPerFace;
  FACE_COLS = config.numColumnsPerFace;
  NUM_SPECIES = config.numSpecies;
  NUM_ESTEPS = config.numEnergySteps;
  NUM_MUSTEPS = config.numMuSteps;
  NUM_OBS = config.numObservers;
  AdiabaticChangeAlg = config.adiabaticChangeAlg;
  AdiabaticFocusAlg = config.adiabaticFocusAlg;

  TOTAL_NUM_SHELLS = config.numNodesPerStream;

  config.simStartTimeDay = config.simStartTime / DAY;
  config.simStopTimeDay  = config.simStopTime / DAY;
  config.tDel            /= DAY;
  config.masInitTimeStep /= DAY;

  config.mhdUs           = ( config.flowMag / C );
  config.mhdNsAu         = ( config.mhdDensityAu / MHD_DENSITY_NORM );
  config.mhdBsAu         = ( config.mhdBAu / MHD_B_NORM );

}

#!/bin/bash

trap ctrl_c INT

function ctrl_c() {
  ${echo} "${cR}==> Caught CTRL-C, shutting down!${cX}"
  exit 1
}

function display_help {
echo "____ ___  ____ ____ _  _    ___ ____ ____ ___    ____ _  _ _ ___ ____ "
echo "|___ |__] |__/ |___ |\/|     |  |___ [__   |     [__  |  | |  |  |___ "
echo "|___ |    |  \ |___ |  |     |  |___ ___]  |     ___] |__| |  |  |___ "
echo "                                                                      "
echo ""
echo "      TEST SUITE v1.0.0"
echo "USAGE:   ./run_test_suite.sh"
echo ""
echo "By default, the above command will run the test suite on all"
echo "default tests using the 'eprem' executable from the '../bin/'"
echo "folder (assuming it has been built)."
echo ""
echo "OPTIONS:"
echo ""
echo "Flag options with choices are specified as '-opt=<OPTION>'."
echo ""
echo "-nochecksetup   Don't check the environment."
echo ""
echo "-epremexe=        Use this to run the testsuite on a specific eprem executable."
echo "                This should be a full path and is useful for development tests."
echo ""
echo "-mpicall=       Use this to run the testsuite with a specific mpi calling mechanism."
echo "                The call should end with '-np' or equivalent as the # of ranks will"
echo "                be placed after the call. "
echo "                The default is 'mpirun -np'."
echo ""
echo "-np=            Number of MPI processes to use.  Note that, due to low resolution,"
echo "                the number of processes is limited.  Default is 1."
echo ""
echo "-test=          Comma-seperated list of subset of tests to run."
echo "                Also can be used to run non-standard/experimental tests."
echo ""
echo "-nocleanup      By default, only the initial and final conditions of a run are "
echo "                kept in the run folders. Set this to keep the full run."
echo ""
echo "-clean          Clear out all testsuite runs."
echo ""
echo "-norun          Does not run the EPREM code.  Checks for a previous run and compares if found."
echo ""
echo "-nocompare      Do not compare the runs to their reference runs."
echo ""
echo "-compareprec=   Set the precision for the comparisons (decimal place)"
echo "                Default is 5."
echo ""
echo "-nocolor        Disable color text output."
echo ""
}

#Set number of processors to use for testsuite:
np=6
norun=0
nocompare=0
nocleanup=0
clean=0
nochecksetup=0
nocolor=0
setrefdata=0
epremexe="eprem"
mpicall="mpirun -np"

AVAIL_TEST_RUNS_LIST="seed
wind
cone
mas_wind
mas_cme"

TEST_RUNS_LIST=${AVAIL_TEST_RUNS_LIST}

for i in "$@"
do
case $i in
    -norun)
    norun=1
    ;;
    -nocompare)
    nocompare=1
    ;;
    -nocleanup)
    nocleanup=1
    ;;
    -nochecksetup)
    nochecksetup=1
    ;;
    -nocolor)
    nocolor=1
    ;;
    -setrefdata)
    setrefdata=1
    ;;
    -test=*)
    TEST_RUNS_LIST="${i#*=}"
    TEST_RUNS_LIST="${TEST_RUNS_LIST//','/' '}"
    ;;
    -epremexe=*)
    epremexe="${i#*=}"
    ;;
    -mpicall=*)
    mpicall="${i#*=}"
    ;;
    -np=*)
    np="${i#*=}"
    ;;    
    -clean)
    clean=1
    norun=1
    nocompare=1
    nochecksetup=1
    ;;
    -h)
    display_help
    exit 0
    ;;
    --help)
    display_help
    exit 0
    ;;    
    *)
    echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
    echo "ERROR!  Unknown option: $i"
    echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
    display_help
    exit 1    
    ;;
esac
done

# ****** Get test suite parameters ******
WD=${PWD}
RESULTSDIR=${WD}/results
BINDIR=${WD}/../bin
SRCDIR=${WD}/../src
ROOTDIR=${WD}/..
TSLOG=${RESULTSDIR}/testsuite.log

if [ ${nocolor} == 0 ]
then
  cX="\033[0m"
  cR="\033[1;31m"
  cB="\033[1;34m"
  cG="\033[32m"
  cC="\033[1;96m"
  cM="\033[35m"
  cY="\033[1;93m"
  Bl="\033[1;5;96m"
  echo="echo -e"
else
  cX=
  cR=
  cB=
  cG=
  cC=
  cM=
  cY=
  Bl=
  echo="echo"
fi

${echo} " "
${echo} "${cB}%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%${cX}"
${echo} "${cY}____ ___  ____ ____ _  _    ___ ____ ____ ___    ____ _  _ _ ___ ____ ${cX}"
${echo} "${cY}|___ |__] |__/ |___ |\/|     |  |___ [__   |     [__  |  | |  |  |___ ${cX}"
${echo} "${cY}|___ |    |  \ |___ |  |     |  |___ ___]  |     ___] |__| |  |  |___ ${cX}"
${echo} "${vY}                                                                      ${cX}"
${echo} "${cB}%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%${cX}"
${echo} "Welcome to the EPREM test suite!"
#
# ****** Test for correct prerequisites and environment ******
#
${echo} "Checking file structure..."

#
# ****** Test for correct prerequisites and environment ******
#
if [ ${nochecksetup} == 0 ]
then
  ${echo} "Checking software requirements..."
  #Check that python is installed:
  PTEST=$(which python3)
  if [ -z "${PTEST}" ]
  then
    ${echo} "${cR}==> ERROR! Python3 does not seem to be installed!${cX}"
    ${echo} "${cR}    This testsuite requires Python3 with the packages:${cX}"
    ${echo} "${cR}      numpy, netCDF4, jdcal${cX}"
    exit 1
  fi
  ${echo} "${cG}==> Python is installed!${cX}"
 #
 # Check for required packages.
 #
  PYTHON_PKG_LIST="numpy
  jdcal
  netCDF4
  "

  for pypkg in $PYTHON_PKG_LIST
  do
    python3 -c "import ${pypkg}" 2>/dev/null
    pychk=$?
    if [ $pychk -eq 1 ]; then
      ${echo} "${cR}==> ERROR! Missing required package ${pypkg}.${cX}"
      ${echo} "${cR}    This testsuite requires Python3 with the packages:${cX}"
      ${echo} "${cR}      numpy, h5py${cX}"      
      exit 1
    fi
  done
  ${echo} "${cG}==> All required python3 packages are present!${cX}"
#
# Check that the eprem bin directory is in the user's path, if not, add it.
#
  ${echo} "==> Checking PATH...."
  PTEST=$(which epremdigest.py)
  if [ -z "${PTEST}" ]
  then
    ${echo} "${cY}==> WARNING: EPREM not in PATH!${cX}"
    ${echo} "${cY}==> Appending ${BINDIR} to PATH...${cX}"
    export PATH="${BINDIR}:${PATH}"
  fi
  PTEST=$(which epremdigest.py)
  if [ -z "${PTEST}" ]; then
    ${echo} "${cR}==> ERROR! EPREM bin PATH problem!${cX}"
    exit 1
  fi
${echo} "${cG}==> Everything seems OK to run the EPREM testsuite!${cX}"
fi
#
# ****** Check that user is running the tests and really wants to set ref data:
#
if [ ${setrefdata} -eq 1 ] && [ ${norun} -eq 1 ]
then
  ${echo} "${cR} ==> ERROR!  You are trying to set reference data without running the tests!${cX}"
  exit 1
fi
if [ ${setrefdata} -eq 1 ]
then
  ${echo} "${cR}╔═╗┌─┐┌┬┐  ╔╗╔┌─┐┬ ┬  ╦═╗┌─┐┌─┐┌─┐┬─┐┌─┐┌┐┌┌─┐┌─┐  ╔╦╗┌─┐┌┬┐┌─┐${cX}"
  ${echo} "${cR}╚═╗├┤  │   ║║║├┤ │││  ╠╦╝├┤ ├┤ ├┤ ├┬┘├┤ ││││  ├┤    ║║├─┤ │ ├─┤${cX}"
  ${echo} "${cR}╚═╝└─┘ ┴   ╝╚╝└─┘└┴┘  ╩╚═└─┘└  └─┘┴└─└─┘┘└┘└─┘└─┘  ═╩╝┴ ┴ ┴ ┴ ┴${cX}"
  read -p "==> Setting new reference data after runs complete...are you SURE? (y/n):" yn
  if [ ${yn} = "y" ]
  then
    read -p "==> Are you really really really SURE?! (y/n):" yn
    if [ ${yn} = "n" ]
    then
      ${echo} "==> ${cR}Aborting!${cX}"
      exit 0
    fi
  else
    ${echo} "==> ${cR}Aborting!${cX}"
    exit 0
  fi
fi

#
# ****** Remove previous comparison results if making new ones.
#
if [ ${nocompare} == 0 ]
then
  rm ${RESULTSDIR}/files/* 2>/dev/null
  rm ${RESULTSDIR}/files/* 2>/dev/null
  rm ${RESULTSDIR}/eprem_testsuite_results.txt 2>/dev/null
fi

########################################################################
########################################################################
##
## ****** Start loop through test problems ******
##
########################################################################
########################################################################

Ti=0
for TESTNAME in ${TEST_RUNS_LIST}
do
  Ti=$((${Ti}+1))
#
# ****** Make sure test is in the testsuite:
#
  ${echo} "${cC}%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%${cX}"
  testok=0
  for run_test in ${AVAIL_TEST_RUNS_LIST}; do
    if [[ ${run_test} == ${TESTNAME} ]]; then
      testok=1
      break
    fi
  done
  if [ ${testok} -eq 0 ]
  then
    ${echo} "${cR}==> TEST ${cX}${cM}${TESTNAME}${cX}${cR} is not a valid test in the testsuite!${cX}"
    exit 1
  fi
  ${echo} "STARTING TEST ${cM}${TESTNAME}${cX}"
  ${echo} "${cC}%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%${cX}"
  ${echo} "==> Gathering test information..."
# ****** Get directories:
  RUNDIR=${WD}/${TESTNAME}/run
  REFDIR=${WD}/${TESTNAME}/reference
  INPUTDIR=${WD}/${TESTNAME}/input

  if [ ! -d ${RUNDIR} ]
  then
    ${echo} "${cR}!!!> ERROR! Run directory does not exist for test ${TESTNAME}!${cX}"
    exit 1
  fi

  if [ ! -d ${REFDIR} ]
  then
    ${echo} "${cR}!!!> ERROR! Reference directory does not exist for test ${TESTNAME}!${cX}"
    exit 1
  fi

  if [ ! -d ${INPUTDIR} ]
  then
    ${echo} "${cR}!!!> ERROR! Input directory does not exist for test ${TESTNAME}!${cX}"
    exit 1
  fi

  if [ ! -e ${INPUTDIR}/eprem_input_file ]
  then
    ${echo} "${cR}!!!> ERROR! Test ${TESTNAME} does not have an input file!${cX}"
    exit 1
  fi

  cd ${RUNDIR}
  
  if [ ${clean} == 1 ]
  then
    ${echo} "==> Clearing run directory..."
    rm -fr ${RUNDIR}/*
  fi

  if [ ${norun} == 0 ]
  then

    if [ -e ${RUNDIR}/epremRunParams.dat ]
    then
      ${echo} "==> Clearing old run..."
      rm -fr ${RUNDIR}/*
    fi

    ${echo} "==> Copying input files..."
    cp -r ${INPUTDIR}/* ${RUNDIR}/ 2>/dev/null

    cd ${RUNDIR}

    ${echo} "======================================================="
    ${echo} "${cB}==> RUNNING EPREM${cX}"
    ${echo} "======================================================="
    ${echo} "==> Running eprem with command:"
    ${echo} "==> ${mpicall} ${np} ${epremexe} eprem_input_file 1>eprem.log 2>eprem.err"
    ${mpicall} ${np} ${epremexe} eprem_input_file 1>eprem.log 2>eprem.err

    # Check that a completed run exists in the run folder
    run_completed_test=$(grep "RUN COMPLETE" ${RUNDIR}/eprem.log)
    if [ -z "${run_completed_test}" ]
    then
      if [ ${norun} == 0 ]
      then
        ${echo} "${cR}!!!> ERROR! Test ${TESTNAME} did not run correctly!${cX}"
        ${echo} "${cR}!!!> Check the run folder: ${RUNDIR} ${cX}"
        ${echo} "${cR}!!!> eprem.log contents: ${cX}"
        cat ${RUNDIR}/eprem.log
        ${echo} "${cR}!!!> eprem.err contents: ${cX}"
        cat ${RUNDIR}/eprem.err
        exit 1
      fi
    fi
  
    if [ ${setrefdata} -eq 1 ] && [ ${norun} -eq 0 ]
    then
      ${echo} "${cR}=======================================================${cX}"
      ${echo} "${cR}==> SETTING REFERENCE DATA FOR RUN${cX}"
      ${echo} "${cR}=======================================================${cX}"

      ${echo} "${cR}==> Removing old reference data...${cX}"
      rm -fr ${REFDIR}/* 2>/dev/null
      ${echo} "${cR}==> Copying current run data into reference directory...${cX}"
      cp ${RUNDIR}/* ${REFDIR}/
      ${echo} "${cR}==> Generating text digest of netCDF output...${cX}"
      epremdigest.py ${RUNDIR} -o ${REFDIR}/eprem_digest.txt
      rm -f ${REFDIR}/*.nc
    fi
#
# ****** Compare run data.
#
    if [ ${nocompare} == 0 ]
    then
      ${echo} "======================================================="
      ${echo} "${cB}==> COMPARING RUN TO REFERENCE DATA${cX}"
      ${echo} "======================================================="
      ${echo} "==> Running comparison..."
      ${echo} "==> Generating digest of run output..."
      result_file=${RESULTSDIR}/files/${TESTNAME}\_eprem\_compare\_run.txt
      rm $result_file 2>/dev/null
      touch $result_file
      
      echo "=======================================================" >> $result_file
      echo "${TESTNAME} RESULTS" >> $result_file
      echo "=======================================================" >> $result_file
      epremdigest.py ${RUNDIR} -o ${RUNDIR}/eprem_digest.txt
      if [ ! -e ${REFDIR}/eprem_digest.txt ]
      then
        echo "ERROR! No reference digest found at ${REFDIR}/eprem_digest.txt" >> $result_file
        echo "FAIL" >> $result_file
      else
        if diff -u ${REFDIR}/eprem_digest.txt ${RUNDIR}/eprem_digest.txt >> $result_file 2>&1
        then
          echo "PASS" >> $result_file
        else
          echo "FAIL" >> $result_file
        fi
      fi
      #See if run passed:
      PASS_FAIL[${Ti}]=$(tail -n 1 $result_file)
      passfailcomp=( ${PASS_FAIL[${Ti}]} )
      if [ "${passfailcomp[0]}" = "FAIL" ]
      then
        ${echo} "${cR}==> Test seems to have FAILED!${cX}"
      else
        ${echo} "${cG}==> Test seems to have PASSED!${cX}"
      fi
      ${echo} "==> Adding comparison to summary file..."
      cat $result_file >> ${RESULTSDIR}/eprem\_testsuite\_results.txt
    fi
#
# ****** Clean run data.
#
    if [ ${nocleanup} == 0 ]
    then
      if [ ${norun} == 0 ]
      then
        ${echo} "======================================================="
        ${echo} "${cB}==> CLEANING RUN DATA${cX}"
        ${echo} "======================================================="
        ${echo} "==> Removing files from run data..."
        rm -fr ${RUNDIR}/*
      fi
    fi
  fi
done
${echo} "${cC}%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%${cX}"


# Display summary and timing results.
if [ ${nocompare} == 0 ]
then
  ${echo} "${cY}==================================================${cX}"
  ${echo} "${cY}Summary of test results:${cX}"
  ${echo} "${cY}==================================================${cX}"
  echo "==================================================" > ${RESULTSDIR}/passfail.txt
  echo "Summary of test results:" >> ${RESULTSDIR}/passfail.txt
  echo "==================================================" >> ${RESULTSDIR}/passfail.txt

  Ti=0
  for TESTNAME in ${TEST_RUNS_LIST}
  do
    Ti=$((${Ti}+1))
    passfailcomp=( ${PASS_FAIL[${Ti}]} )
    if [ "${passfailcomp[0]}" = "FAIL" ]; then
      ${echo} "$(printf "%-25s %-25s" ${TESTNAME} ${cR}"${PASS_FAIL[${Ti}]}"${cX})"
      echo "$(printf "%-25s %-25s" ${TESTNAME} "${PASS_FAIL[${Ti}]}")">> ${RESULTSDIR}/passfail.txt
    else
      ${echo} "$(printf "%-25s %-25s" ${TESTNAME} ${cG}"${PASS_FAIL[${Ti}]}"${cX})"
      echo "$(printf "%-25s %-25s" ${TESTNAME} "${PASS_FAIL[${Ti}]}")" >> ${RESULTSDIR}/passfail.txt
    fi
  done
  echo "=======================================================" >> ${RESULTSDIR}/passfail.txt
  cat ${RESULTSDIR}/eprem_testsuite_results.txt >> ${RESULTSDIR}/passfail.txt
  mv ${RESULTSDIR}/passfail.txt ${RESULTSDIR}/eprem_testsuite_results.txt
fi

${echo} "${cC}%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%${cX}"

exit 0

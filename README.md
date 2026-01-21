# EPREM-PSI: A version of EPREM 
    
[Predictive Science Inc.](https://www.predsci.com)  
 
--------------------------------  
  
## OVERVIEW ##
EPREM-PSI is...  
  
 
--------------------------------  
   
## HOW TO BUILD EPREM-PSI ##
  
EPREM-PSI has been tested to work using GCC's `gfortran`, Intel's `ifx`, and NVIDIA's `nvfortran` compilers.
It is recommended to use the latest compiler version available.  

EPREM-PSI requires the [HDF5](https://www.hdfgroup.org/solutions/hdf5), NETCDF, and LIBCONFIG libraries.  
  
1. Find the closest configuration file in `/conf` and make a local copy of it.  
2. Modify the conf file to set the library paths/flags and compiler flags to match your system environment.  
3. Run the build script (for example, `./build.sh <your_custom_config.conf>`).  
4. It is recommended to add the `bin` folder to your system path.  
  
### RUN THE EPREM-PSI TESTSUITE ###
  
To test if the installation is working, we recommend running the testsuite after installation.  
To do this, enter the `testsuite/` directory and run:  
  
`./run_test_suite.sh`  
  
This will run the tests with `bin/eprem` using 6 MPI ranks.  
  
NOTE:  Currently the testsuite has no validation tests.
  
--------------------------------  
  
## HOW TO RUN EPREM-PSI ##
  
### Setting Input Options  
  
`EPREM-PSI` uses a text input file.  The default name for the input is `eprem_input_file`  
  
  
### Launching the Code ###
    
To run `EPREM-PSI`, set the desired run parameters in the input file, then copy or link the `eprem` executable into the same directory as the input file and run the command:  
  
`<MPI_LAUNCHER> <MPI_OPTIONS> -np <N> ./eprem <input_file>`  
  
where `<N>` is the total number of MPI ranks to use, `<MPI_LAUNCHER>` is your MPI run command (e.g. `mpiexec`,`mpirun`, `ibrun`, `srun`, etc), and `<MPI_OPTIONS>` are additional MPI options that may be needed (such as `--bind-to socket` or `--bind-to numa` for CPUs running with OpenMPI).  

For example:  `mpirun -np 1 ./eprem eprem_input_file`  
  
  
### Solution Output ###
  
The output of EPREM-PSI are netcdf files.
  
### Processing Scripts ###
  
The `/bin` folder contains python and bash scripts that can be used with EPREM-PSI.  

--------------------------------

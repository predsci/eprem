# EPREM-PSI: A version of EPREM 
    
[Predictive Science Inc.](https://www.predsci.com)  
 
--------------------------------  
  
## OVERVIEW ##
EPREM-PSI is...  
  
 
--------------------------------  
   
## HOW TO BUILD EPREM-PSI ##
  
EPREM-PSI has been tested to work using GCC's `gfortran` (>8), Intel's `ifx`, and NVIDIA's `nvfortran` compilers.
It is recommended to use the latest compiler version available.  

EPREM-PSI requires the HDF4, NETCDF, LIBCONFIG, and [HDF5](https://www.hdfgroup.org/solutions/hdf5) libraries.  
The libraries should be compiled by the same compiler HipFT is using (e.g. nvfortran).  
  
1. Find the build script from the `build_examples` folder that is closest to your setup and copy it into the top-level directory.  
2. Modify the script to set the library paths/flags and compiler flags compatible with your system environment.  
3. Modify the script to set the compiler options to reflect your setup.  
4. Run the build script (for example, `./my_build.sh`).  
5. It is recommended to add the `bin` folder to your system path.  
  
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

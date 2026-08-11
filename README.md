# EPREM
    
[Predictive Science Inc.](https://www.predsci.com)  
 
--------------------------------  
  
## OVERVIEW ##
`EPREM` is ...  
This version includes the capability to couple with [`MAS`](https://github.com/predsci/mas) simulation output.  
It is part of `STAT`.  
It can run on `CORHEL-CME` simulation output.  
 
--------------------------------  
   
## HOW TO BUILD EPREM ##
  
`EPREM` has been tested to work using GCC's `gfortran`, Intel's `ifx`, and NVIDIA's `nvfortran` compilers.
It is recommended to use the latest compiler version available.  

`EPREM` requires the [HDF5](https://www.hdfgroup.org/solutions/hdf5), NETCDF, and LIBCONFIG libraries.  
  
1. Find the closest configuration file in `/conf` and make a local copy of it.  
2. Modify the conf file to set the library paths/flags and compiler flags to match your system environment.  
3. Run the build script (for example, `./build.sh <your_custom_config.conf>`).  
4. It is recommended to add the `bin` folder to your system path.  
  
### RUN THE EPREM TESTSUITE ###
  
To test if the installation is working, we recommend running the testsuite after installation.  
To do this, enter the `testsuite/` directory and run:  
  
`./run_test_suite.sh`  
  
This will run the tests with `bin/eprem` using 6 MPI ranks.  
  
NOTE:  Currently the testsuite does not validate the solutions, it just runs the tests.
  
--------------------------------  
  
## HOW TO RUN EPREM ##
  
### Setting Input Options  
  
`EPREM` uses a text input file.  The default name for the input is `eprem_input_file`  
  
  
### Launching the Code ###
    
To run `EPREM`, set the desired run parameters in the input file, then copy or link the `eprem` executable into the same directory as the input file and run the command:  
  
`<MPI_LAUNCHER> <MPI_OPTIONS> -np <N> ./eprem <input_file>`  
  
where `<N>` is the total number of MPI ranks to use, `<MPI_LAUNCHER>` is your MPI run command (e.g. `mpiexec`,`mpirun`, `ibrun`, `srun`, etc), and `<MPI_OPTIONS>` are additional MPI options that may be needed (such as `--bind-to socket` or `--bind-to numa` for CPUs running with OpenMPI).  

For example:  `mpirun -np 1 ./eprem eprem_input_file`  
  
  
### Solution Output ###
  
The output of EPREM are netcdf files.
  
### Processing Scripts ###
  
The `/bin` folder contains python and bash scripts that can be used with EPREM.  

--------------------------------

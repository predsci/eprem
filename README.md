<img height=100 src="doc/eprem_logo.png" alt="EPREM"/>  

# Energetic Particle Radiation Environment Model #
   
[Predictive Science Inc.](https://www.predsci.com)  
 
--------------------------------  
  
## OVERVIEW ##
The Energetic Particle Radiation Environment Model (`EPREM`) is a physics-based solar energetic particle (SEP) transport code.   It models energetic particle acceleration and transport
using a Lagrangian system, which comoves with the plasma.  Details of the model can be found [here](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2009SW000523).  
This version of EPREM includes the capability to couple with [`MAS`](https://github.com/predsci/mas) MHD simulation output.  
EPREM is a core component of the SPE Threat Assessment Tool ([`STAT`](https://iopscience.iop.org/article/10.1088/1742-6596/1225/1/012007)).  
It can run on [`CORHEL-CME`](https://ccmc.gsfc.nasa.gov/models/CORHEL-CME~1) MAS simulation output.  
 

### Please Read First: ###
EPREM is a complex code with many years of development and contains numerous input parameters, use cases, experimental features, deprecated features, etc.   The EPREM code provided here will become part of future STAT open-source release as a sub-module.  Until then, it is primarily intended for reference purposes and/or for advanced users/developers.  
  
==> Please contact PSI at support@predsci.com if you plan to use this code for research purposes. <==

--------------------------------  
   
## HOW TO INSTALL EPREM ##
  
### Compilers ###
EPREM has been tested to work with the following compilers:  
  
 - GCC's `gfortran`  
 - NVIDIA's `nvfortran`   
 - INTEL's `ifx`  

It is recommended to use the latest compiler version available.  

### Dependencies ###

`EPREM` requires the [HDF5](https://www.hdfgroup.org/solutions/hdf5), [NETCDF](https://www.unidata.ucar.edu/software/netcdf), and [LIBCONFIG](https://hyperrealm.github.io/libconfig) libraries.  
  
### Build Instructions ###

1. Find the configuration file in `/conf` that is closest to your system aand make a local copy of it.  
2. Modify the copied conf file to set the library paths/flags and compiler flags to match your system environment.  
3. Run the build script (for example, `./build.sh <your_custom_config.conf>`).  
4. It is recommended to add the `bin` folder to your system path.  
  
### Run the Testsuite ###
  
To ensure the installation was successful, a testsuite is provided in the `testsuite` folder.  
To run the testsuite:
```
> cd testsuite
> ./run_test_suite.sh
```
This will run the tests with `bin/eprem` using 6 MPI ranks.  
The available options can be viewed by running `./run_test_suite.sh -h`.  
  
NOTE:  Currently the testsuite does not validate the solutions, it just runs the tests.
  
--------------------------------  
  
## HOW TO RUN EPREM ##
  
### Setting Input Options ###
  
`EPREM` uses a text input file.  The default name for the input is `eprem_input_file`  
Examples of input file scan be found in the testsuite's "input" folders.
  
  
### Launching the Code ###
    
To run `EPREM`, set the desired run parameters in the input file and run the command (assuming eprem is in your path):  
  
`<MPI_LAUNCHER> <MPI_OPTIONS> -np <N> eprem <input_file>`  
  
where `<N>` is the total number of MPI ranks to use, `<MPI_LAUNCHER>` is your MPI run command (e.g. `mpiexec`,`mpirun`, `ibrun`, `srun`, etc), and `<MPI_OPTIONS>` are additional MPI options that may be needed (such as `--bind-to socket` or `--bind-to numa` for CPUs running with OpenMPI).  

For example:  `mpirun -np 6 eprem eprem_input_file`  
  
  
### Solution Output ###
  
The output of EPREM are netcdf files, one for each stream that contains either the distribution function or differential fluxes (depending on options) for all nodes on the stream for all time outputs.
  
### Processing Scripts ###
  
The `bin` folder contains python and bash scripts that can be used to read, analyze, and plot eprem simulation output.  

--------------------------------

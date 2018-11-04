# Code experiments
This folder contains all applications mentioned in the bachelor thesis "Using the Legion Programming Framework to extract Implicit Parallelism with Logical Regions"

## Dependencies
* [PAPI](http://icl.cs.utk.edu/projects/papi/wiki/Introduction_to_PAPI-C)
* working Legion/Regent Installation
* Python 2.X Installation to use Legion Prof and Legion Spy

## 0_SHP_Dokulil (not published)
This folder contains the SPH-Simulation, that was provided by Dr. Jiri Dokulil. The original execution time measurements were replaced with PAPI.

Build with `make` and execute with `make run`. Building with `make sequential` results in a version that does not make use of OpenMP

## Regent Applications

All regent applications share a common makefile, launch parameters and possibilities to replace the default values in the code.

The paths to `regent.py`, `legion_prof.py` and `legion_spy.py` must be specified in the makefile of each regent application. Also the number of CPUs to launch the application on, as this variable is needed by the Legion runtime.

Run with `make run`.  
Run Legion Prof with `make prof`  
Run Legion Spy with `make spy`  
Or all at once with `make all`  
  
`make clean` removes the logfiles used by Legion Prof and Legion Spy.  
`make clean-all` removes the logfiles and all generated profiling and Spy output!  

Profiling data is saved in the folder `./profiler_html/<DATETIME_of_Run>`  
Legion Spy data is saved in the folder
`./profiler_html/spy_graphs/<DATETIME_of_Run>`  
  
launching a regent application with the `-h` flag prints all possible launch parameters:
  
  * `-c` -- set particle count
  * `-i` -- set iteration count
  * `-m` -- set the space boundaries to [-max;max]

Additionally every regent application reads its default variables from the file `sph_config.rg`, where those hardcoded values can easily be changed permanently, when using the launch parameters seems unpractical.

### 1_SPH_Regent_Naive
contains a naive translation of the original SPH simulation with every C++ array as its own region. The execution of the program is sequential as it does not use partitioning.

### 2_SPH_Regent_SingleR
contains the program with all regions consolidated into one single region. The execution of the program is still sequential.

### 3_SPH_Regent_Partitioning
contains the first version of the program to use partitioning. However it produces incorrect results as this version does not declare correct region requirements. Comparing the dependency graphs of this implementation with Regent_implicit gives a good example of how to debug such errors.

### 4_SPH_Regent_implicit
contains a correct implementation of the simulation, that can be executed in parallel/distributed.

### 5_SPH_Regent_explicitGhost
contains a correct implementation of the simulation, that can be executed in parallel/distributed. But introduces a full copy of all particles for ghost access with an explicit synchronization scheme. While only Regent operations are used to achieve this, explicit synchronization breaks with Legion's programming model.

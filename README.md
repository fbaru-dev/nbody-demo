# Demo Session for Intel® Advisor and Intel® Compiler C++
This is an example code based on a simple N-body simulation of a distribution of point masses placed
at location r_1,...,r_N and have masses m_1,...,m_N. The position of the particles after a specified
time is computed using a finite difference methods for ordinary differential equation.

## Implementation
For each particle the position, the velocity, the acceleration and the mass is stored in a C-like
structure and for an N particles case, an array of this structure is allocated. This is the 
simple data-structure which is very close to the physical representation of a particle mass.
The file `Particle.hpp` contains the implementation of such data-structure.

For each particle indexed by i, the accelearation is computed a_i = G*mj*(ri-rj)/|ri-rj|^3, which 
value is used to update the velocity and position using the Euler integration scheme.
Furthermore the total energy of the particles' group is computed.
The file `GSimulation.cpp` contains the implementation of the algorithm.

## Directory structure of the Demo
The demo consists of several directories, which correspond to the different
optimization steps to take to enabling vectorization and OpenMP multi-threding of the code.
Each directory has its onw makefile to compile and run the test case.
To compiler the code type `make` and the run the simulation type `make run`.
As benchmark, the simulation starts with 2000 particles and 500 integration steps. One can
change the default giving the number of particles and the number of integration steps using
the command line argument:
`./nbody.x < # of particles> < # of integration>`

Try to change the number of particles and observe how the performance changes.

## Different versions
To start the demo, go to the folder `ver0`, compile and run the test.

### Intial version: ver0
The typical output of the simulation is:
```
Run the default test case on CPU: 
 ./nbody.x 
 ===============================
  Initialize Gravity Simulation
  nPart = 2000; nSteps = 500; dt = 0.1
 ------------------------------------------------
  s       dt      kenergy     time (s)    GFlops      
 ------------------------------------------------
   50      5      103.29      1.5109      3.84        
  100     10      389.23      1.3833      4.1941      
  150     15      461.36      1.3327      4.3534      
  200     20      902.23      1.4129      4.1064      
  250     25      959.14      1.4233      4.0764      
  300     30      1042.5      1.5151      3.8293      
  350     35      1438        1.4837      3.9103      
  400     40      1567.5      1.306       4.4424      
  450     45      2000        1.419       4.0888      
  500     50      2666.5      1.3435      4.3184         

 # Number Threads     : 1
 # Total Time (s)     : 14.466
 # Average Perfomance : 4.0513 +- 0.11807
```

On output is printed some useful information. Colomnwise: s is the
number of steps; dt is the physical time taking into account the physical
time integration step; kenery is the kinetic energy of the group of particles;
time is the computational time taken till that time step; GFlops is the
number of giga flops per second. 
N.B. The GFlops is an estimation done by looking into the code and counting
the number of math operations according to the algorithm. This is used only
as standard metric for comparison. More realistic numbers can be measured
in different way (Roofline model of Intel® Advisor).

Following the five steps of code modernization, 
https://software.intel.com/en-us/articles/what-is-code-modernization
we can improve the performance of the code.

- describe the Intel® Advisor result
- compile the code with processor specific optimization: -xSSE4.2, -xAVX, -xCORE-AVX2, -xCORE-AVX512, -xMIC-AVX512
- generate the compiler report and describe the different options: 
-  -qopt-report[=N]: default level is 2
-  -qopt-report-phase=<vec,loop,openmp,...>: default is all
-  -qopt-report-file=stdout | stderr | filename
-  -qopt-report-filter="GSimulation.cpp,130-204"

Then show how verbose is the compiler report and use filtering.

### ver1 
Solution of the ver0. The optimization are: -O2 -xAVX or higher. 
The Makefile is the only difference. Here we generate higher vectorized code and
produce the compiler report.
One should run this version in the same way as before and:
- show the new performance numbers
- describe the Intel® Advisor result
- generate the compiler report
- explain FP conversions and precision of constants, variables and math functions

### ver2
Solution of the ver1. The difference is in the GSimulation.cpp file where the consistent
computation with floats is made (constants and SQRT function).
One should run this version in the same way as before and:
- show the new performance numbers
- describe the Intel® Advisor result
- generate the compiler report
- explain the remark #25085: Preprocess Loopnests: Moving Out Load and Store and 
  remark #15415: vectorization support: non-unit strided load was generated for the variable
  ....
  remark #15300: LOOP WAS VECTORIZED
  remark #15452: unmasked strided loads: 6 
  remark #15475: --- begin vector cost summary ---
  remark #15476: scalar cost: 115 
  remark #15477: vector cost: 26.750 
  remark #15478: estimated potential speedup: 4.070 
  remark #15488: --- end vector cost summary ---
  ....
- explain vectorization gather/scatter
- explain AoS and SoA differences

### ver3
Solution of the ver2. The differences are in:
- Particle.hpp: the new SoA data structure is implemented
- GSimulation.hpp: modified the data member according to SoA
- GSimulation.cpp: allocation and reference to SoA

One should run this version in the same way as before and:
- show the new performance numbers
- describe the Intel® Advisor result
- generate the compiler report
- explain the remark #15344: loop was not vectorized: vector dependence prevents vectorization
  remark #15346: vector dependence: assumed ANTI dependence between ... and ...
  remark #15346: vector dependence: assumed FLOW dependence between ... and ...
- explain the vectorization and how much we gain using it
- refer to the Intel® compiler autovectorization guide and explain the requirements
  for autovectorization
- explain #pragma simd
- explain #pragma simd reduction
- modify in the `Makefile` the CXXFLAGS adding the OMPFLAGS at line 8, recompile and run
  remark #15301: OpenMP SIMD LOOP WAS VECTORIZED
- at this point running the code shows wrong results (Warning with SIMD, be aware of the full control)
- try to use #pragma simd reduction (solution in the file GSimulation-simd.cpp)
  NB rember that the simd reduction is not allowed on `particles->acc_x[i]`
  Solution: 
    - cp GSimulation.cpp GSimulation.cpp.bkp
    - cp GSimulation-simd.cpp GSimulation.cpp
  recompile and run
  remark #15301: OpenMP SIMD LOOP WAS VECTORIZED
- rerun and show that the result is now correct

### ver4
This is the clean solution of the ver3 after all modification done live in the
previous session.
One should run this version in the same way as before and:
- show the new performance numbers
- describe the Intel® Advisor result
- generate the compiler report
- explain remark #15389: vectorization support: reference ... has `unaligned` access 
- explain the data alignment with examples and the alignment size (16/32/64 bytes)
- exlpain peel and reminder loops

## ver5
This is the solution of the ver4, with all the allocations replaced by the memory
alignment allocation function.
Running this version allows to see that even modifing the memory allocation functions,
the data is not aligned. One needs to use the function `__assume_aligned(...)`.
Recompile the code adding the option: -DASALIGN.
One should run again this version with the alignment option and:
- show the new performance numbers
- describe the Intel® Advisor result
- generate the compiler report

This concludes the basic vectorization part of the demo.
At this point, only two topics are missing:
- advanced cache optimization (loop-tiling) (ver6)
- enabling OpenMP (ver7)

### ver6
This is the cache optimized version of the code, without OpenMP.
The performance depends on the size of the tile and the number or particles.
One should run again this version and:
- describe in detail what is this kind of optimization and how depends on the tile size
- show the new performance numbers
- describe the Intel® Advisor result
- generate the compiler report

### ver7
This is the version of the code with OpenMP. Play with the number of threads,
openmp scheduling and threads affinity.

### ver8
This is the version of the code with OpenMP and cache tiling.
One can also play with the floating point model -fp-model fast=2, for example and
look for further performance improvements

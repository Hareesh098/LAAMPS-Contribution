Dear Harish,

As promised, I start sharing my code to study oscillatory instability. My plan is to share all the 
codes in three steps.
Step 1: codes to perform athermal quasi static simulations AQS.
Step2: codes to compute Jacobian
Step3: codes to diagonalize Jacobian

In this email, I have attached all the codes for Step 1 with a brief explanation. You might have many 
questions and doubts, feel free to email me.

Thanks,
Joyjit

------------------------------------------------------------------------------------------------------
# generating initial configuration
------------------------------------------------------------------------------------------------------
We start simulations from a configuration where a binary assembly randomly arranged in a fixed square 
box.

1.Compile: g++ iconf-granfrictionBinary.cc -o iconf-granfrictionBinary
2.Run: ./iconf-granfrictionBinary 500 0.93 2 0.5 0.7 1

500: nb_particles
0.93: packing fraction
2: dimension
0.5: radius of small particle
0.7: radius of large particle
1: ratio between small particles and large particles (1.0 for binary mixture)

The above code will generate an initial configuration inside file initconf.d
-------------------------------------------------------------------------------------------------------

-------------------------------------------------------------------------------------------------------
# including new force field in LAMMPS
-------------------------------------------------------------------------------------------------------
1.Download LAMMPS version 11Aug17 (The latest version should also work but it was not tested so far).
2.Go to the directory: lammps-11Aug17/src
3.Compile LAMMPS: make serial or make mpi

If installation is successful: 
4. Copy the following two files in lammps-11Aug17/src :
    pair_gran_hertzpolynomial_history.h and pair_gran_hertzpolynomial_history.cpp
5. Go to the directory: lammps-11Aug17/src
6. make clean-all
7. make serial or make mpi

At this end, LAMMPS successfully incorporates the new force field and we are ready to run numerical 
simulations.
-------------------------------------------------------------------------------------------------------

-------------------------------------------------------------------------------------------------------
# performing simulations using LAMMPS
-------------------------------------------------------------------------------------------------------
To perform a typical athermal quasi static simulation AQS:
1. Copy in.aqs and initconf.d in the working directory
2. Run: ~/lammps-11Aug17/src/lmp_serial -in in.aqs -log log.aqs
or
mpirun -np 4 ~/lammps-11Aug17/src/lmp_mpi -in in.aqs -log log.aqs

In a single processor this will take typically ~2mins to finish the job.

If the code runs successfully, it will produce:
1. two dump files: dump-mu-10.lammpstrj and dumplocal-mu-10.lammpstrj.
First one contains particle positions at strain value 0 (i.e. just before AQS). Second one contains 
information about pairs (different components of force and shear).

2.a file called look0.001 which contains information about:
# 1. step 2. strain 3. temperature 4. pressure 5. shear stress_xx 6. _yy 7._xy

3. Two dumpfiles: dump-s-0.001-ss-10.lammpstrj and dumplocal-s-0.001-ss-10.lammpstrj. These are 
information about particle positions and pairs in the sheared frame when strain value is 0.001.
-------------------------------------------------------------------------------------------------------

-------------------------------------------------------------------------------------------------------
# how to read output files
-------------------------------------------------------------------------------------------------------
Any dumplocal-* file contains information about pairs in this order:
1. pair-index, 
2. atom1 id of the pair, 
3. atom2 id of the pair, 
4. x component of normal force, 
5. y component of normal force, 
6. ratio between tangential distance and tangential threshold, 
7. x component of tangential force, 
8. y component of tangential force, 
9. tangential displacement along x, 
10. tangential displacement along y, 
11. x component of pair vector r_{ij} i.e. dx
12. y component of pair vector r_{ij} i.e. dy

Extracting various information associated with tangential force is very important to compute Jacobian. 
If you look at Eqs. 12 and 13 in the supplemental of our recently submitted paper (attached), you 
probably better understand the columns 6, 7, 8, 9 and 10.

These are standard and can be found in LAMMPS manual. ix and iy are integer values stand for 
image_x and image_y of that particle. These terms help us to determine a particle which was initially 
inside the box, where it displaces at time t in the Euclidean space (not necessarily inside the 
simulation box). For orthogonal box [Lx, Ly], the relation would be:
x_true = x + ix * Lx
y_true = y + iy * Ly

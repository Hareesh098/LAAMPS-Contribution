**10June 2025**
========================================================================================================================
Added the following four file based on gran_hooke and gran_hooke_history:
- pair_gran_harmonic.cpp
- pair_gran_harmonic.h
- pair_gran_harmonic_history.cpp
- pair_gran_harmonic_history.h

All it does is following two things:
- Normalise the normal Hookean normal force by sum of two discs interacting:
ccel = ccel =  kn * ((radsum-r)/radsum) - DampingNormal
- And remove the meff from the damping forces: DampingNormal and DampingTangential. So,  in the above,
in normal direction: 
DampingNormal =   damp = gamman * vnnr * rsqinv
and in the tangential dirction:
DampingTangential = gammat * vtr1
[Compiled with]
make yes-molecule yes-bpm yes-extra-molecule yes-extra-pair yes-granular yes-misc yes-brownian yes-colloid yes-dipole
yes-extra-fix yes-extra-compute yes-mc yes-opt
make -j 32 serial
make -j 32 mpi
========================================================================================================================

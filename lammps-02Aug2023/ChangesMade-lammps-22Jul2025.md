**02 August 2023 Stable 4**
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

[Configured with]
make yes-molecule yes-asphere yes-bpm yes-brownian yes-cg-spica yes-class2 yes-colloid yes-dipole yes-dpd-basic 
yes-dpd-meso yes-extra-compute yes-extra-dump yes-extra-fix yes-extra-molecule yes-extra-pair yes-granular yes-mc 
yes-misc yes-ptm yes-reaction yes-rigid yes-shock yes-tally yes-uef

[Compiled with]
make -j 32 serial
make -j 32 mpi
========================================================================================================================

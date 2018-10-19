## Eoin Murray
## 110367833

This program simulates annealing of a noble gas fcc atomic cluster with a monte-carlo molecular-dynamics approach.
The atoms are confined by the Leonard-Jones potential and are annealed by adding energy and integrating the equations of motion.
The heat capacity of the system is investigated by adding incremental energy to the system and allowing the atoms to gain a brownian velocity.
Further systematics of the vibrational spectrum are investigated.

The code runs through the sections of the report available at [http://eoinmurray.io/molec].

Code compiles with

	gfortran dran1_.o ffc_posn.c randveladd.c forces.c int_dyn.c _fft1.o code.c o- binary

None but code.c are provided here.

#!/bin/bash
## template for running MPI programs with qsub for SGE1
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
## The following can also be passed on the qsub command line, perhaps
## more convenient for multiple runs with different numbers of PEs.
#$ -pe orte 4
#
MPI_DIR=/opt/openmpi

$MPI_DIR/bin/mpirun -np $NSLOTS some-program >some-outputfile

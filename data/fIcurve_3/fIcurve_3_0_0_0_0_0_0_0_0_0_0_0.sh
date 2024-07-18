#!/bin/bash
#$ -cwd
#$ -N my_batch
#$ -q cpu.q
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -l h_rt=00:30:00
#$ -o /home/andres/PVcellSimulation/data/fIcurve_3/fIcurve_3_0_0_0_0_0_0_0_0_0_0_0.run
#$ -e /home/andres/PVcellSimulation/data/fIcurve_3/fIcurve_3_0_0_0_0_0_0_0_0_0_0_0.err

source ~/.bashrc
mpiexec -n $NSLOTS -hosts $(hostname) nrniv -python -mpi sim/init.py simConfig=data/fIcurve_3/fIcurve_3_0_0_0_0_0_0_0_0_0_0_0_cfg.json netParams=data/fIcurve_3/fIcurve_3_netParams.py

        
"""
init.py

Starting script to run NetPyNE-based PT model.

Usage:
    python sim/init.py # Run simulation, optionally plot a raster

MPI usage:
    mpiexec -n 4 nrniv -python -mpi sim/init.py
"""

#import matplotlib; matplotlib.use('Agg')  # to avoid graphics error in servers

from netpyne import sim

#cfg, netParams = sim.loadFromIndexFile('index.npjson')
"""
init.py

Starting script to run NetPyNE-based M1 model.

Usage:
    python init.py # Run simulation, optionally plot a raster

MPI usage:
    mpiexec -n 4 nrniv -python -mpi init.py

Contributors: salvadordura@gmail.com
"""

import matplotlib; matplotlib.use('Agg')  # to avoid graphics error in servers
from netpyne import sim


#------------------------------------------------------------------------------
## Function to modify cell params during sim (e.g. modify PT ih)


# -----------------------------------------------------------
# Main code
# Option to run one example
#cfg, netParams = sim.loadFromIndexFile('index.npjson')
# Option necessary for batch to work
cfg, netParams = sim.readCmdLineArgs(simConfigDefault='src/cfg.py', netParamsDefault='src/netParams.py')

sim.createSimulateAnalyze(netParams, cfg)

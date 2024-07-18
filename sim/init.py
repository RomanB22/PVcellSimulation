"""
init.py

Starting script to run NetPyNE-based PT model.

Usage:
    python sim/init.py # Run simulation, optionally plot a raster

MPI usage:
    mpiexec -n 4 nrniv -python -mpi sim/init.py
"""

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
import matplotlib.pyplot as plt
import numpy as np

#------------------------------------------------------------------------------
## Function to modify cell params during sim (e.g. modify PT ih)
# -----------------------------------------------------------
# Main code
# Option to run one example
#cfg, netParams = sim.loadFromIndexFile('index.npjson')
# Option necessary for batch to work
cfg, netParams = sim.readCmdLineArgs(simConfigDefault='sim/cfg.py', netParamsDefault='sim/netParams.py')

sim.createSimulateAnalyze(netParams, cfg)

# Access the recorded spike data
spike_data = sim.allSimData['spkt'] 
print(len(spike_data))

# Parameters
time_window = [0, 2000]  # Define the time window for the PSTH
bin_size = 10  # Define the bin size for the PSTH
Min_trials=50


record_pops = [f'PV5B_{i}' for i in range(1, 50)]
spike_times, spike_ids = sim.allSimData['spkt'], sim.allSimData['spkid']  # Spike times and corresponding cell IDs

# Filter spike data for populations that start with 'PV5B_'
pv5b_spike_times = []
cell_ids = [cell.gid for cell in sim.net.cells if cell.tags['pop'] in record_pops]

for (j,spike_id) in enumerate(spike_ids):
    if spike_id in cell_ids:
        pv5b_spike_times.append(spike_times[j])
    
# pv5b_spike_times = spike_times[np.isin(spike_ids, cell_ids)]


bins = np.arange(time_window[0], time_window[1] + bin_size, bin_size)
spike_counts = np.zeros(len(bins) - 1)

# Bin the spike times
for neuron_spike_times in pv5b_spike_times:
    spike_counts += np.histogram(neuron_spike_times, bins=bins)[0]

# Normalize by the number of neurons and bin size to get spike density
spike_density = spike_counts / (Min_trials * bin_size)

# Plot the TDFR
plt.figure(figsize=(10, 5))
plt.bar(bins[:-1], spike_density, width=bin_size, align='edge', color='black')
plt.xlabel('Time (s)')
plt.ylabel('Spike Density (spikes/s)')
plt.title('Time Dependent Firing Rate')
plt.savefig('tdfr.png') 




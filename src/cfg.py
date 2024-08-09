"""
cfg.py

Simulation configuration for M1 model (using NetPyNE)
"""

from netpyne import specs
import pickle
import numpy as np

cfg = specs.SimConfig()

#------------------------------------------------------------------------------
#
# SIMULATION CONFIGURATION
#
#------------------------------------------------------------------------------

# current injection params
ExpAmps = list(np.arange(0.08, 0.54, 0.02))  # amplitudes
ExpTargetRates = [1.14, 3., 4.85, 6.57, 7.71, 9.14, 10.43, 11.43, 12.43, 12.86, 13.86, 14.43, 15.15, 15.86, 16.29,
				  17.29, 17.57, 18.57, 19.14, 19.57, 19.43, 19.71, 20.2, 20.]
amps = ExpAmps[:8]  # [::3]
amps.insert(0, 0.04)
amps.insert(0, 0.02)
times = list(np.arange(1000, 2000 * len(amps), 2000))  # start times
dur = 400  # ms
targetRates = ExpTargetRates[:8]  # [::3]
targetRates.insert(0, 0)
targetRates.insert(0, 0)
targetRates = [i * 1000 / 400 for i in targetRates]


#------------------------------------------------------------------------------
# Run parameters
#------------------------------------------------------------------------------

cfg.duration = 2000 * len(amps)
cfg.dt = 0.1
cfg.seeds = {'conn': 4321, 'stim': 1234, 'loc': 4321}
cfg.hParams = {'celsius': 34, 'v_init': -80}
cfg.verbose = False
cfg.createNEURONObj = True
cfg.createPyStruct = True
cfg.cvode_active = False
cfg.cvode_atol = 1e-6
cfg.cache_efficient = True
cfg.printRunTime = 0.1
cfg.includeParamsLabel = True
cfg.printPopAvgRates = True
cfg.checkErrors = True
cfg.connRandomSecFromList = False

#------------------------------------------------------------------------------
# Recording
#------------------------------------------------------------------------------
# record from all cells
cfg.recordTraces = {'V_soma': {'sec':'soma', 'loc':0.5, 'var':'v'}}
cfg.recordSpikesGids = ['all']
cfg.recordStim = True
cfg.recordTime = True
cfg.recordStep = 0.1
cfg.recordLFP = False #[[10, y, 90] for y in range(450, 1250, 100)]
cfg.saveLFPCells = False

#------------------------------------------------------------------------------
# Saving
#------------------------------------------------------------------------------
cfg.saveFolder = 'data/'
cfg.savePickle = False
cfg.saveJson = True
cfg.saveDataInclude = ['simData', 'simConfig', 'netParams', 'net']
cfg.backupCfgFile = None #['cfg.py', 'backupcfg/']
cfg.gatherOnlySimData = False
cfg.saveCellSecs = False
cfg.saveCellConns = True

#------------------------------------------------------------------------------
# Cells
#------------------------------------------------------------------------------
cfg.Experiment = 'fI' # 'fI' or 'PT5B_inputs', to distinguish whether we run sims to calibrate in-vitro fI curve or to reproduce in-vivo inputs

#------------------------------------------------------------------------------
# Current inputs
#------------------------------------------------------------------------------
if cfg.Experiment == 'fI':
	cfg.addIClamp = True
	cfg.addVecStim = False
	cfg.addNetStim = False
	# current injection params
	cfg.IClamp1 = {'pop': 'FoxP2', 'sec': 'soma', 'loc': 0.5, 'dur': dur, 'amp': amps, 'start': times}
	cfg.simLabel = 'FoxP2_fI/FoxP2'
#------------------------------------------------------------------------------
# VecStim inputs
#------------------------------------------------------------------------------
cfg.Go = 'Go'
cfg.Condition = 'InVivo' + '_%s' % cfg.Go  # 'InVivo', 'OnlyIncre', 'MirrorDecre'
cfg.ESynMech = 'AMPA'  # ['AMPA', 'NMDA']
cfg.AMPANMDAWeightsIncre = 0.005
cfg.AMPANMDAWeightsDecre = 0.005
cfg.AMPANMDAWeightsNotChanging = 0.005  # or 0.005/2 depending on if we want 400 pA of AMPA current or 200 pA per presynaptic spike
cfg.delay = 1.0
cfg.preStim = 1800
cfg.postStim = 1800
cfg.Random = 0
cfg.RangeConnections = [37, 90] #[55, 87] # Minimum and maximum number of connections from CFA PT5B to FoxP2
# Not sure yet how to cut the upper limit. Not used yet
cfg.FoxP2 = True # Whether to use the calibrated LTS3.hoc model or to use the FS3.hoc

#------------------------------------------------------------------------------
# Synapses
#------------------------------------------------------------------------------
cfg.AMPATau2Factor = 5.0 # In order to match an EPSC peak of 400 pA. Divide by 2 if you want to have 200 pA as maximum current
cfg.somaProb = 0.2 # Probability of connection to the soma. Extracted from Ach inputs from motor neurons to FoxP2

#------------------------------------------------------------------------------
# Analysis and plotting
#------------------------------------------------------------------------------
with open('cells/popColors.pkl', 'rb') as fileObj: popColors = pickle.load(fileObj)['popColors']

if cfg.Experiment == 'fI':
	cfg.analysis['plotfI'] = {'amps': amps, 'times': times, 'dur': dur, 'target': {'rates': targetRates}, 'saveFig': True, 'showFig': False, 'calculateFeatures': ''}
	cfg.analysis['plotTraces'] = {'include': ['FoxP2'], 'timeRange': [0,cfg.duration], 'oneFigPer': 'cell', 'figSize': (10,4), 'saveFig': True, 'showFig': False}

if cfg.Experiment == 'PT5B_inputs':
	#####################
	cfg.addNetStim = True # Add the rest of physiological inputs to FoxP2
	cfg.NetStimRate = 0.0001 # From firing rate in Hz to Interval the conversion is Interval[ms] = 1000/Freq[Hz]
	cfg.NetStimNoise = 0.5 # Fraction of noise in NetStim (0 = deterministic; 1 = completely random)
	cfg.NetStimWeight = 0.005
	cfg.NetStimNumber = 1e10 # Max number of spikes generated (default = 1e12)
	cfg.NetStimDelay = 1
	cfg.simLabel = 'FoxP2_VecStim_%s/FoxP2_' % cfg.Condition if cfg.FoxP2 else 'PV_VecStim_%s/PV_' % cfg.Condition
	cfg.duration = cfg.preStim + cfg.postStim
	#####################
	cfg.addIClamp = True # I clamp to simulate the change in resting potential in-vivo
	cfg.IAmp = 0 # nA
	# current injection params
	cfg.IClamp1 = {'pop': 'FoxP2', 'sec': 'soma', 'loc': 0.5, 'dur': cfg.duration, 'amp': cfg.IAmp, 'start': 0}
	#####################
	cfg.addVecStim = True
	#####################
	timeRange = [200, cfg.duration-200]
	cfg.analysis['plotTraces'] = {'include': [('FoxP2', i) for i in range(5)], 'timeRange': timeRange,
								  'oneFigPer': 'trace', 'overlay': False, 'figSize': (10, 15), 'saveFig': True,
								  'showFig': False}
	cfg.analysis['plotRaster'] = {'include': ['FoxP2'], 'timeRange': timeRange, 'orderInverse': False, 'saveFig': True,
								  'showFig': False}
	cfg.analysis['plotSpikeFreq'] = {'include': ['FoxP2'], 'timeRange': timeRange, 'measure': 'rate', 'binSize': 25,
									 'saveFig': True, 'showFig': False, 'density' : False,
									 'xlabel' : 'Time (ms)', 'marker' : 'x'}

#------------------------------------------------------------------------------
# Parameters
#------------------------------------------------------------------------------
# example of how to set params; but set from batch.py
cfg.tune = specs.Dict()
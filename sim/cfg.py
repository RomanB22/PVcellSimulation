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

#------------------------------------------------------------------------------
# Run parameters
#------------------------------------------------------------------------------
cfg.duration = 1*1e3
cfg.dt = 0.05
cfg.seeds = {'conn': 4321, 'stim': 1234, 'loc': 4321}
cfg.hParams = {'celsius': 34, 'v_init': -80}
cfg.verbose = True
cfg.createNEURONObj = True
cfg.createPyStruct = True
cfg.cvode_active = False
cfg.cvode_atol = 1e-6
cfg.cache_efficient = True
cfg.printRunTime = 0.1

cfg.includeParamsLabel = False
cfg.printPopAvgRates = False

cfg.checkErrors = True

#------------------------------------------------------------------------------
# Recording
#------------------------------------------------------------------------------
allpops = ['PV5B']

cfg.recordTraces = {'V_soma': {'sec':'soma', 'loc':0.5, 'var':'v'}}

cfg.recordStim = True
cfg.recordTime = True
cfg.recordStep = 0.1
cfg.recordLFP = False #[[10, y, 90] for y in range(450, 1250, 100)]
cfg.saveLFPCells = False

#------------------------------------------------------------------------------
# Saving
#------------------------------------------------------------------------------
cfg.simLabel = 'PVcell'
cfg.saveFolder = '.'
cfg.savePickle = False
cfg.saveJson = True
cfg.saveDataInclude = ['simData', 'simConfig', 'netParams']#, 'net']
cfg.backupCfgFile = None #['cfg.py', 'backupcfg/']
cfg.gatherOnlySimData = False
cfg.saveCellSecs = False
cfg.saveCellConns = False

#------------------------------------------------------------------------------
# Cells
#------------------------------------------------------------------------------
cfg.ihModel = 'migliore'  # ih model
cfg.ihGbar = 1.0  # multiplicative factor for ih gbar in PT cells
cfg.ihGbarBasal = 1.0 # 0.1 # multiplicative factor for ih gbar in PT cells
cfg.ihlkc = 0.2 # ih leak param (used in Migliore)
cfg.ihlkcBasal = 1.0
cfg.ihlkcBelowSoma = 0.01
cfg.ihlke = -86  # ih leak param (used in Migliore)
cfg.ihSlope = 14*2

cfg.removeNa = False  # simulate TTX; set gnabar=0s
cfg.somaNa = 5
cfg.dendNa = 0.3
cfg.axonNa = 7
cfg.axonRa = 0.005

cfg.gpas = 0.5  # multiplicative factor for pas g in PT cells
cfg.epas = 0.9  # multiplicative factor for pas e in PT cells

#------------------------------------------------------------------------------
# Synapses
#------------------------------------------------------------------------------
cfg.synWeightFractionEE = [0.5, 0.5] # E->E AMPA to NMDA ratio
cfg.synWeightFractionEI = [0.5, 0.5] # E->I AMPA to NMDA ratio
cfg.synWeightFractionSOME = [0.9, 0.1] # SOM -> E GABAASlow to GABAB ratio

cfg.AMPATau2Factor = 1.0

#------------------------------------------------------------------------------
# Network
#------------------------------------------------------------------------------
cfg.weightNormThreshold = 4.0  # weight normalization factor threshold

#------------------------------------------------------------------------------
# Subcellular distribution
#------------------------------------------------------------------------------
cfg.addSubConn = False

#------------------------------------------------------------------------------
# Current inputs
#------------------------------------------------------------------------------
cfg.addIClamp = True

# current injection params
amps = list(np.arange(0.0, 0.65, 0.05))  # amplitudes
times = list(np.arange(1000, 2000 * len(amps), 2000))  # start times
dur = 500  # ms
targetRates = [0., 0., 19., 29., 37., 45., 51., 57., 63., 68., 73., 77., 81.]

cfg.IClamp1 = {'pop': 'PV5B', 'sec': 'soma', 'loc': 0.5, 'dur': dur, 'amp': amps, 'start': times}


#------------------------------------------------------------------------------
# NetStim inputs
#------------------------------------------------------------------------------
cfg.addNetStim = False

cfg.NetStim1 = {'pop': 'PV5B', 'ynorm':[0,1], 'sec': 'apic_5', 'loc': 0.5, 'synMech': ['AMPA'], 'synMechWeightFactor': [1.0],
				'start': 0, 'interval': 1000.0/40.0, 'noise': 0.0, 'number': 1000.0, 'weight': 10.0, 'delay': 0}

#------------------------------------------------------------------------------
# Analysis and plotting
#------------------------------------------------------------------------------
with open('cells/popColors.pkl', 'rb') as fileObj: popColors = pickle.load(fileObj)['popColors']

cfg.analysis['plotfI'] = {'amps': amps, 'times': times, 'dur': dur, 'target': {'rates': targetRates}, 'saveFig': True, 'showFig': True, 'calculateFeatures': ''}
cfg.analysis['plotTraces'] = {'include': [('PV5B',0)], 'timeRange': [0,cfg.duration], 'oneFigPer': 'cell', 'figSize': (10,4), 'saveFig': True, 'showFig': False}
#cfg.analysis['plotLFP'] = {'separation': 1.0, 'plots': ['timeSeries', 'locations'], 'saveFig': True, 'showFig': False}

#------------------------------------------------------------------------------
# Parameters
#------------------------------------------------------------------------------
# example of how to set params; but set from batch.py
cfg.tune = specs.Dict()
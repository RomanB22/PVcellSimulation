"""
batch.py 

Batch simulation for M1 model using NetPyNE

Contributors: salvadordura@gmail.com
"""
import numpy as np
from netpyne import specs
from netpyne.batch import Batch

""" Not used here, but may be useful if we modify the model
# ----------------------------------------------------------------------------------------------
# Weight Normalization Exc
# ----------------------------------------------------------------------------------------------

#pops =  ['IT2', 'IT4', 'IT5A', 'IT5B', 'PT5B', 'IT6', 'CT6', 'PV2', 'SOM2'],
def weightNormE(pops = ['PV5B'], secs = None, locs = None,
                allSegs = True, rule = 'PV_reduced', weights =list(np.arange(0.01, 0.2, 0.01)/100.0)):

    # Add params
    from netParams import netParams
    from cfg import cfg

    excludeSegs = ['axon']
    if not secs:
        secs = []
        locs = []
        for secName, sec in netParams.cellParams[rule]['secs'].items():
            if secName not in excludeSegs:
                if allSegs:
                    nseg = sec['geom']['nseg']
                    for iseg in list(range(nseg)):
                        secs.append(secName) 
                        locs.append((iseg+1)*(1.0/(nseg+1)))
                else:
                    secs.append(secName) 
                    locs.append(0.5)

    params = specs.ODict()
    params[('NetStim1', 'pop')] = pops
    params[('NetStim1', 'loc')] = locs
    params[('NetStim1', 'sec')] = secs
    params[('NetStim1', 'weight')] = weights

    groupedParams = [('NetStim1', 'sec'), ('NetStim1', 'loc')]


    initCfg = {}
    initCfg['duration'] = 1.0*1e3
    initCfg[('analysis','plotTraces','timeRange')] = [0, 1000]
    initCfg['weightNorm'] = False
    initCfg['stimSubConn'] = False
    initCfg['addNetStim'] = True
    initCfg[('NetStim1', 'synMech')] = ['AMPA','NMDA']
    initCfg[('NetStim1','synMechWeightFactor')] = [0.5,0.5]
    initCfg[('NetStim1', 'start')] = 700
    initCfg[('NetStim1', 'interval')] = 1000
    initCfg[('NetStim1','ynorm')] = [0.0, 1.0]

    initCfg[('NetStim1', 'pop')] = ['PV5B']

    initCfg[('NetStim1', 'noise')] = 0
    initCfg[('NetStim1', 'number')] = 1
    initCfg[('NetStim1', 'delay')] = 1
    #initCfg[('GroupNetStimW1', 'pop')] = 'None'
    initCfg['addIClamp'] = 0

    b = Batch(params=params, netParamsFile='sim/netParams.py', cfgFile='sim/cfg.py', initCfg=initCfg, groupedParams=groupedParams)

    return b


# ----------------------------------------------------------------------------------------------
# EPSPs via NetStim
# ----------------------------------------------------------------------------------------------
def EPSPs():
    params = specs.ODict()

    params['groupWeight'] = [x*0.05 for x in np.arange(1, 8, 1)]
    params['ihGbar'] = [0.0, 1.0]
 
    
    # initial config
    initCfg = {}
    initCfg['duration'] = 0.5*1e3
    initCfg['addIClamp'] = False
    initCfg['addNetStim'] = True
    initCfg[('GroupNetStimW1', 'pop')] = 'PV5B'
    initCfg[('analysis','plotTraces','timeRange')] = [0, 500]
    initCfg['excTau2Factor'] = 2.0
    initCfg['weightNorm'] = True
    initCfg['stimSubConn'] = False
    initCfg['ihGbarZD'] = None

    groupedParams = [] 

    b = Batch(params=params, netParamsFile='sim/netParams.py', cfgFile='sim/cfg.py', initCfg=initCfg,
              groupedParams=groupedParams)

    return b
"""

# ---------------------------------------------------------------------------------------------- #
# -----                              f-I curve calibration                                ------ #
# ---------------------------------------------------------------------------------------------- #
def evolCellPV5B():
    # parameters space to explore
    params = specs.ODict()

    TotalCellSurface = 3969.4023 # um2
    ExpCapacitance = 47.7 # pF
    SpecificCapacitance = ExpCapacitance*1e-6/(TotalCellSurface*1e-8) # uF/cm2
    ExpInputResistance = 239 # MOhm
    SpecificLeakConductance = 1/(239*1e6)/(TotalCellSurface*1e-8) # S/cm2

    params[('tune', 'soma', 'Ra')] = [150.*0.5, 150*1.5]
    params[('tune', 'soma', 'cm')] = [SpecificCapacitance*0.5, SpecificCapacitance*1.5]
    params[('tune', 'soma', 'kdrin', 'gkdrbar')] = [0.018*0.5, 0.018*1.5]
    params[('tune', 'soma', 'hin', 'gbar')] = [0.00001*0.5, 0.00001*0.5]
    params[('tune', 'soma', 'kapin', 'gkabar')] = [0.0032*15*0.5, 0.0032*15*0.5]
    params[('tune', 'soma', 'Nafx', 'gnafbar')] = [0.045*0.5, 0.045*1.5]
    params[('tune', 'soma', 'pas', 'e')] = [-73*1.5, -73.0*0.5]
    params[('tune', 'soma', 'pas', 'g')] = [SpecificLeakConductance*0.5, SpecificLeakConductance*1.5]
    params[('tune', 'soma', 'L')] = [27 * 0.5, 27 * 1.5]

    params[('tune', 'dend', 'Ra')] = [150.*0.5, 150*1.5]
    params[('tune', 'dend', 'cm')] = [SpecificCapacitance*0.5, SpecificCapacitance*1.5]
    params[('tune', 'dend', 'kdrin', 'gkdrbar')] = [0.018*0.5*0.5, 0.018*0.5*1.5]
    params[('tune', 'dend', 'kapin', 'gkabar')] = [0.0032*15*10*0.5, 0.0032*15*10*0.5]
    params[('tune', 'dend', 'Nafx', 'gnafbar')] = [0.018*5*0.5, 0.018*5*1.5]
    params[('tune', 'dend', 'pas', 'e')] = [-73*1.5, -73.0*0.5]
    params[('tune', 'dend', 'pas', 'g')] = [SpecificLeakConductance*0.5, SpecificLeakConductance*1.5]

    # current injection params
    ExpAmps = list(np.arange(0.08, 0.54, 0.02))  # amplitudes
    ExpTargetRates = [1.14, 3., 4.85, 6.57, 7.71, 9.14, 10.43, 11.43, 12.43, 12.86, 13.86, 14.43, 15.15, 15.86, 16.29, 17.29, 17.57, 18.57,$
    amps = ExpAmps[::3]
    amps.insert(0, 0.04)
    amps.insert(0, 0.02)
    times = list(np.arange(1000, 2000 * len(amps), 2000))  # start times
    dur = 500  # ms
    targetRates = ExpTargetRates[::3]#[0., 0., 19., 29., 37., 45., 51., 57., 63., 68., 73., 77., 81.]
    targetRates.insert(0,0)
    targetRates.insert(0,0)
    # initial cfg set up
    initCfg = {} # specs.ODict()
    initCfg['duration'] = 2000 * len(amps)
    initCfg[('hParams', 'celsius')] = 37

    initCfg['savePickle'] = False
    initCfg['saveJson'] = True
    initCfg['saveDataInclude'] = ['simConfig', 'netParams', 'net', 'simData']

    initCfg[('IClamp1', 'pop')] = 'PV5B'
    initCfg[('IClamp1', 'amp')] = amps
    initCfg[('IClamp1', 'start')] = times
    initCfg[('IClamp1', 'dur')] = 1000
    initCfg[('analysis', 'plotTraces')] = {'include': [('PV5B', 0)], 'timeRange': [0, initCfg['duration']],
                                           'oneFigPer': 'cell', 'figSize': (10, 4),
                                           'saveFig': True, 'showFig': False}

    initCfg[('analysis', 'plotfI', 'amps')] = amps
    initCfg[('analysis', 'plotfI', 'times')] = times
    initCfg[('analysis', 'plotfI', 'dur')] = dur
    initCfg[('analysis', 'plotfI', 'target')] = {'rates': targetRates}

    for k, v in params.items():
        initCfg[k] = v[0]  # initialize params in cfg so they can be modified

    # fitness function
    fitnessFuncArgs = {}
    fitnessFuncArgs['target'] = {'rates': targetRates}

    def fitnessFunc(simData, **kwargs):
        targetRates = kwargs['target']['rates']
        diffRates = [abs(x-t) for x,t in zip(simData['fI'], targetRates)]
        fitness = np.mean(diffRates)

        print(' Candidate rates: ', simData['fI'])
        print(' Target rates:    ', targetRates)
        print(' Difference:      ', diffRates)

        return fitness


    # create Batch object with paramaters to modify, and specifying files to use
    b = Batch(cfgFile='sim/cfg.py', netParamsFile='sim/netParams.py', params=params, initCfg=initCfg)

    # Set evol method (all param combinations)
    b.method = 'evol'
    b.evolCfg = {
        'evolAlgorithm': 'custom',
        'fitnessFunc': fitnessFunc, # fitness expression (should read simData)
        'fitnessFuncArgs': fitnessFuncArgs,
        'pop_size': 4,
        'num_elites': 1, # keep this number of parents for next generation if they are fitter than children
        'mutation_rate': 0.4,
        'crossover': 0.5,
        'maximize': False, # maximize fitness function?
        'max_generations': 30,
        'time_sleep': 50, # wait this time before checking again if sim is completed (for each generation)
        'maxiter_wait': 20, # max number of times to check if sim is completed (for each generation)
        'defaultFitness': 1000 # set fitness value in case simulation time is over
    }

    return b
# ----------------------------------------------------------------------------------------------
# f-I curve
# ----------------------------------------------------------------------------------------------
def fIcurve(): #TODO: Change values for the IClamp1 so it respects the experimental protocol
    params = specs.ODict()

    params[('IClamp1', 'pop')] = ['PV5B']
    params[('IClamp1', 'amp')] = list(np.arange(0.0, 10.0, 1)/10.0)

    # initial config
    initCfg = {}
    initCfg['duration'] = 1.0*1e3
    initCfg['addIClamp'] = True
    initCfg['addNetStim'] = False
    initCfg['weightNorm'] = True
    initCfg[('IClamp1','sec')] = 'soma'
    initCfg[('IClamp1','loc')] = 0.5
    initCfg[('IClamp1','start')] = 200
    initCfg[('IClamp1','dur')] = 600
    initCfg[('analysis','plotTraces','timeRange')] = [0, 1000]
    initCfg[('analysis', 'plotfI')] = {} # Don't plot fI in this case

    groupedParams = []

    b = Batch(params=params, netParamsFile='sim/netParams.py', cfgFile='sim/cfg.py', initCfg=initCfg, groupedParams=groupedParams)
    b.method = 'grid'

    return b

# ----------------------------------------------------------------------------------------------
# Run configurations
# ----------------------------------------------------------------------------------------------
def setRunCfg(b, type='mpi_bulletin', nodes=1, coresPerNode=8):
    if type=='mpi_bulletin':
        b.runCfg = {'type': 'mpi_bulletin',
            'script': 'sim/init.py',
            'skip': False,
            'skipCustom': '_data.json'}

    elif type=='mpi_direct':
        b.runCfg = {'type': 'mpi_direct',
            'cores': 1,
            'script': 'sim/init.py',
            'mpiCommand': 'mpiexec',
            'skip': False,
            'skipCustom': '_data.pkl'}

# ----------------------------------------------------------------------------------------------
# Main code
# ----------------------------------------------------------------------------------------------

if __name__ == '__main__': 

    #b = fIcurve()
    #b.batchLabel = 'fIcurve'
    b = evolCellPV5B()
    b.batchLabel = 'evolfI'
    b.saveFolder = 'data/'+b.batchLabel
    setRunCfg(b, 'mpi_direct')
    b.run() # run batch 

"""
batch.py 

Batch simulation for M1 model using NetPyNE

Contributors: salvadordura@gmail.com
"""
import numpy as np
from netpyne import specs
from netpyne.batch import Batch
import os

""" Not used here, but may be useful if we modify the model
# ----------------------------------------------------------------------------------------------
# Weight Normalization Exc
# ----------------------------------------------------------------------------------------------

#pops =  ['IT2', 'IT4', 'IT5A', 'IT5B', 'PT5B', 'IT6', 'CT6', 'PV2', 'SOM2'],
def weightNormE(pops = ['FoxP2'], secs = None, locs = None,
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

    initCfg[('NetStim1', 'pop')] = ['FoxP2']

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
    initCfg[('GroupNetStimW1', 'pop')] = 'FoxP2'
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


def batchRun(Condition='InVivo', Go='Go', RandomSpikeInit=0, preStim=1800, postStim=1800):
    # parameters space to explore
    params = specs.ODict()
    params[('tune', 'Random')] = RandomSpikeInit
    initCfg = {} # specs.ODict()
    initCfg['preStim'] = preStim
    initCfg['postStim'] = postStim
    initCfg['dt'] = 0.1
    initCfg['recordStep'] = 0.1
    initCfg[('hParams', 'celsius')] = 37
    initCfg['Go'] = Go
    initCfg['Condition'] = Condition
    initCfg['addVecStim'] = True
    initCfg['AMPANMDAWeightsIncre'] = 0.005
    initCfg['AMPANMDAWeightsDecre'] = 0.005
    initCfg['AMPANMDAWeightsNotChanging'] = 0.005

    for k, v in params.items():
        initCfg[k] = v[0]  # initialize params in cfg so they can be modified

    # create Batch object with paramaters to modify, and specifying files to use
    b = Batch(cfgFile='src/cfg.py', netParamsFile='src/netParams.py', params=params, initCfg=initCfg)

    # Set evol method (all param combinations)
    b.method = 'grid'
    return b

# ---------------------------------------------------------------------------------------------- #
# -----                              f-I curve calibration                                ------ #
# ---------------------------------------------------------------------------------------------- #
def evolCellFoxP2():
    # parameters space to explore
    params = specs.ODict()

    TotalCellSurface = 9803.4718 #3969.4023 # um2
    ExpCapacitance = 47.7 # pF
    SpecificCapacitance = 1.2 #ExpCapacitance*1e-6/(TotalCellSurface*1e-8) # uF/cm2
    ExpInputResistance = 239 # MOhm
    SpecificLeakConductance = 6e-5 #1/(239*1e6)/(TotalCellSurface*1e-8) # S/cm2
    #print(SpecificLeakConductance, SpecificCapacitance)
    #[222.15039017222813, 1.4144239759846486, 0.07787265804398133, 0.0, 0.0, 
    #0.6309169430849156, 0.0013728455162340843, -97.45589500657026, 4.890896421208355e-05] 
    """
    [226.8778026684533, 1.4668784547221765, 0.0742666956942066, 8.181467787372347e-07, 9.065582224086196e-06, 0.6393160722589732, 
     0.0013812904703444037, -98.9999992775243, 4.650114503867427e-05]
    """
    params[('tune', 'soma', 'Ra')] = [222*0.95, 222*1.05]
    params[('tune', 'soma', 'cm')] = [1.4*0.95, 1.4*1.05]
    params[('tune', 'soma', 'kdrin', 'gkdrbar')] = [0.077872658*0.95, 0.077872658*1.05]
    params[('tune', 'soma', 'hin', 'gbar')] = [0.000002*0.01, 0.000002*1]
    params[('tune', 'soma', 'kapcb', 'gkabar')] = [0.00075*0.01, 0.0075*1]
    params[('tune', 'soma', 'Nafx', 'gnafbar')] = [0.63091*0.95, 0.63091*1.05]
    params[('tune', 'soma', 'catcb', 'gcatbar')] = [0.0013728*0.95, 0.0013728*1.05]
    params[('tune', 'soma', 'pas', 'e')] = [-99,-40]
    params[('tune', 'soma', 'pas', 'g')] = [4.89089e-5*0.95, 4.89089e-5*1.05]
    #params[('tune', 'soma', 'L')] = [10, 60]


    #params[('tune', 'dend', 'Ra')] = [150.*0.5, 150*1.5]
    #params[('tune', 'dend', 'cm')] = [SpecificCapacitance*0.8, SpecificCapacitance*1.2]
    #params[('tune', 'dend', 'kdrin', 'gkdrbar')] = [0.018*5*0.5, 0.018*5*1.5]
    #params[('tune', 'dend', 'kapcb', 'gkabar')] = [0.00875*0.5, 0.00875*1.5]
    #params[('tune', 'dend', 'Nafx', 'gnafbar')] = [0.018*0.5, 0.018*1.5]
    #params[('tune', 'dend', 'pas', 'e')] = [-64*1.5, -64*0.5]
    #params[('tune', 'dend', 'pas', 'g')] = [SpecificLeakConductance*0.8, SpecificLeakConductance*1.2]

    # current injection params
    ExpAmps = list(np.arange(0.08, 0.54, 0.02))  # amplitudes
    ExpTargetRates = [1.14, 3., 4.85, 6.57, 7.71, 9.14, 10.43, 11.43, 12.43, 12.86, 13.86, 14.43, 15.15, 15.86, 16.29, 17.29, 17.57, 18.57, 19.14, 19.57, 19.43, 19.71, 20.2, 20.]
    amps = ExpAmps[:8] #[::3]
    amps.insert(0, 0.04)
    amps.insert(0, 0.02)
    times = list(np.arange(1000, 2000 * len(amps), 2000))  # start times
    dur = 400  # ms
    targetRates = ExpTargetRates[:8] #[::3]
    targetRates.insert(0,0)
    targetRates.insert(0,0)
    targetRates = [i*1000/400 for i in targetRates]
    #print(amps, targetRates)
    # initial cfg set up
    initCfg = {} # specs.ODict()
    initCfg['duration'] = 2000 * len(amps)
    initCfg[('hParams', 'celsius')] = 37

    initCfg['savePickle'] = False
    initCfg['saveJson'] = True
    initCfg['saveDataInclude'] = ['simConfig', 'netParams', 'net', 'simData']

    initCfg[('IClamp1', 'pop')] = 'FoxP2'
    initCfg[('IClamp1', 'amp')] = amps
    initCfg[('IClamp1', 'start')] = times
    initCfg[('IClamp1', 'dur')] = 400
    initCfg[('analysis', 'plotTraces')] = {'include': [('FoxP2', 0)], 'timeRange': [0, initCfg['duration']],
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
        #if np.mean(simData['fI'])<5: fitness=1000
        #if np.count_nonzero(simData['fI'])<6: fitness=1000
        u, c = np.unique(simData['fI'], return_counts=True)
        dup = u[c > 1]
        if np.count_nonzero(simData['fI'])<4 or np.max(simData['fI'])>100 or sorted(simData['fI']) != simData['fI'] or len(dup)>1: fitness=1000



        print(' Candidate rates: ', simData['fI'])
        print(' Target rates:    ', targetRates)
        print(' Difference:      ', diffRates)

        return fitness


    # create Batch object with paramaters to modify, and specifying files to use
    b = Batch(cfgFile='sim/cfg.py', netParamsFile='sim/netParams.py', params=params, initCfg=initCfg)

    # Set evol method (all param combinations)
    b.method = 'evol'
    b.evolCfg = {
        'evolAlgorithm': 'custom', #'custom',
        'fitnessFunc': fitnessFunc, # fitness expression (should read simData)
        'fitnessFuncArgs': fitnessFuncArgs,
        'pop_size': 100,
        'num_elites': 2, # keep this number of parents for next generation if they are fitter than children
        'mutation_rate': 0.5,
        'crossover': 0.5,
        'maximize': False, # maximize fitness function?
        'max_generations': 100,
        'time_sleep': 500, # wait this time before checking again if sim is completed (for each generation)
        'maxiter_wait': 5, # max number of times to check if sim is completed (for each generation)
        'defaultFitness': 1000 # set fitness value in case simulation time is over
    }

    return b
# ----------------------------------------------------------------------------------------------
# f-I curve
# ----------------------------------------------------------------------------------------------
def fIcurve(): 
    params = specs.ODict()

    # current injection params
    ExpAmps = list(np.arange(0.08, 0.54, 0.02))  # amplitudes
    ExpTargetRates = [1.14, 3., 4.85, 6.57, 7.71, 9.14, 10.43, 11.43, 12.43, 12.86, 13.86, 14.43, 15.15, 15.86, 16.29, 17.29, 
17.57, 18.57, 19.14, 19.57, 19.43, 19.71, 20.2, 20.]
    amps = ExpAmps[:8] #[::3]
    amps.insert(0, 0.04)
    amps.insert(0, 0.02)
    #times = list(np.arange(1000, 2000 * len(amps), 2000))  # start times
    dur = 400  # ms

    Ra, cm, gkdrbar, gbar, gkabar, gnafbar, gcatbar, e, g = 226.8778027, \
    1.4668784547221765,	0.074266696, 8.181467787372347e-07, 9.065582224086196e-06, \
    0.6393160722589732,	0.00138129, -98.99999928, 4.65E-05

    params[('tune','soma', 'Ra')] = [Ra, Ra]
    params[('tune','soma', 'cm')] = [cm, cm]
    params[('tune','soma', 'kdrin', 'gkdrbar')] = [gkdrbar, gkdrbar]
    params[('tune','soma', 'hin', 'gbar')] = [gbar, gbar]
    params[('tune','soma', 'kapcb', 'gkabar')] = [gkabar, gkabar]
    params[('tune','soma', 'Nafx', 'gnafbar')] = [gnafbar, gnafbar]
    params[('tune','soma', 'catcb', 'gcatbar')] = [gcatbar, gcatbar]
    params[('tune','soma', 'pas', 'e')] = [e, e]
    params[('tune','soma', 'pas', 'g')] = [g, g]

    params[('IClamp1', 'pop')] = ['FoxP2']
    params[('IClamp1', 'amp')] = amps
    # initial config
    initCfg = {}
    initCfg['duration'] = 1.0*1e3
    initCfg['addIClamp'] = True
    initCfg['addNetStim'] = False
    initCfg['weightNorm'] = True
    initCfg[('IClamp1','sec')] = 'soma'
    initCfg[('IClamp1','loc')] = 0.5
    initCfg[('IClamp1','start')] = 300
    initCfg[('IClamp1','dur')] = dur

    #initCfg[('analysis','plotTraces','timeRange')] = [0, 1000]
    #initCfg[('analysis', 'plotfI')] = {} # Don't plot fI in this case

    initCfg[('analysis', 'plotTraces')] = {'include': [('FoxP2', 0)], 'timeRange': [0, initCfg['duration']],
                                           'oneFigPer': 'cell', 'figSize': (10, 4),
                                           'saveFig': True, 'showFig': False}

    initCfg[('analysis', 'plotfI', 'amps')] = amps
    initCfg[('analysis', 'plotfI', 'times')] = 300
    initCfg[('analysis', 'plotfI', 'dur')] = dur


    groupedParams = [('tune', 'soma', 'Ra'), ('IClamp1', 'pop'), \
('tune', 'soma', 'cm'), ('tune', 'soma', 'kdrin', 'gkdrbar'), \
('tune', 'soma', 'hin', 'gbar'), ('tune', 'soma', 'kapcb', 'gkabar'), \
('tune', 'soma', 'Nafx', 'gnafbar'), ('tune', 'soma', 'catcb', 'gcatbar'), \
('tune', 'soma', 'pas', 'e'), ('tune', 'soma', 'pas', 'g')]

    b = Batch(params=params, netParamsFile='sim/netParams.py', cfgFile='sim/cfg.py', initCfg=initCfg, groupedParams=groupedParams)
    b.method = 'grid'

    return b

# ----------------------------------------------------------------------------------------------
# Run configurations
# ----------------------------------------------------------------------------------------------
def setRunCfg(b, type='mpi_bulletin', nodes=1, coresPerNode=8):
    if type=='mpi_bulletin':
        b.runCfg = {'type': 'mpi_bulletin',
            'script': 'src/init.py',
            'skip': False,
                    }#'skipCustom': '_data.json'}

    elif type=='mpi_direct':
        b.runCfg = {'type': 'mpi_direct',
            'cores': 1,
            'script': 'src/init.py',
            'mpiCommand': 'mpiexec',
            'skip': False,
            'skipCustom': '_data.pkl'}

    elif type=='hpc_sge':
        b.runCfg = {'type': 'hpc_sge',
                    'jobName': 'my_batch',
                    'cores': 1,
                    'mpiCommand': 'mpiexec',
                    'vmem': '1G',
                    'walltime': "00:30:00",
                    'script': 'src/init.py',
                    'queueName': 'cpu.q',
                    'skip': False,
                    'skipCustom': '_data.json'}


# ----------------------------------------------------------------------------------------------
# Main code
# ----------------------------------------------------------------------------------------------

if __name__ == '__main__':
    #b = fIcurve()
    #b.batchLabel = 'fIcurve'
    #b = evolCellFoxP2()
    #b.batchLabel = 'evolfI_4'
    for Condition in ['InVivo_NoGo','OnlyIncre_NoGo','MirrorDecre_NoGo'
        ,'InVivo_Go','OnlyIncre_Go','MirrorDecre_Go']:
        trialStep = [i for i in range(100)]
        Go = Condition.split('_')[1]
        b = batchRun(Condition=Condition, Go=Go, RandomSpikeInit=[0])
        b.batchLabel = 'Sims_%s' % Condition
        b.saveFolder = 'data/'+b.batchLabel
        setRunCfg(b, 'mpi_bulletin')
        b.run() # run batch

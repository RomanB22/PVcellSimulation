
"""
netParams.py

High-level specifications for M1 network model using NetPyNE
"""
import os
from netpyne import specs
import pickle
import random
import numpy as np

netParams = specs.NetParams()   # object of class NetParams to store the network parameters

netParams.version = 1

try:
    from __main__ import cfg  # import SimConfig object with params from parent module
except:
    from cfg import cfg

#------------------------------------------------------------------------------
#
# NETWORK PARAMETERS
#
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# General connectivity parameters
#------------------------------------------------------------------------------
netParams.defaultThreshold = 0.0 # spike threshold, 10 mV is NetCon default, lower it for all cells
netParams.defaultDelay = 2.0 # default conn delay (ms)
netParams.propVelocity = 500.0 # propagation velocity (um/ms)
netParams.probLambda = 100.0  # length constant (lambda) for connection probability decay (um)
netParams.sizeX = 200
netParams.sizeY = 200  # cortical depth (will be converted to negative values)
netParams.sizeZ = 200

#------------------------------------------------------------------------------
# Cell parameters
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Specification of cell rules not previously loaded
# Includes importing from hoc template or python class, and setting additional params

#------------------------------------------------------------------------------
## PV cell params (3-comp)
if cfg.FoxP2:
    cellRule = netParams.importCellParams(label='FoxP2_reduced', conds={'cellType':'FoxP2', 'cellModel':'HH_reduced'},
    fileName='cells/LTS3.hoc', cellName='FoxP2', cellInstance=True)
    cellRule['secLists']['spiny'] = ['soma', 'dend']

    Ra, cm, e, g, gkdrbar, gbar, gkabar, gnafbar, gcatbar = {}, {}, {}, {}, {}, {}, {}, {}, {}

    Ra['soma'], cm['soma'], e['soma'], g['soma'], gkdrbar['soma'], gbar['soma'], gkabar['soma'], gnafbar['soma'], gcatbar['soma']  = \
    178, 1.46, -81.5, 3.84E-05, 7.13E-02, 8.43E-07, 9.52E-06, 6.11E-01, 1.41E-03

    Ra['dend'], cm['dend'], e['dend'], g['dend'], gkdrbar['dend'], gbar['dend'], gkabar['dend'], gnafbar['dend'], gcatbar['dend'] = \
    178, 1.45, -84, 3.94E-05, 1.91E-01, None, 5.60E-06, 9.59E-01, None

    Ra['axon'], cm['axon'], e['axon'], g['axon'], gkdrbar['axon'], gbar['axon'], gkabar['axon'], gnafbar['axon'], gcatbar['axon'] = \
    172, 1.33, -97, 3.99E-05, 5.36E-01, None, None, 8.77E-01, None

    for sec in cellRule['secs'].keys():
        netParams.cellParams['FoxP2_reduced']['secs'][sec]['geom']['Ra'] = Ra[sec]
        netParams.cellParams['FoxP2_reduced']['secs'][sec]['geom']['cm'] = cm[sec]
        netParams.cellParams['FoxP2_reduced']['secs'][sec]['mechs']['pas']['e'] = e[sec]
        netParams.cellParams['FoxP2_reduced']['secs'][sec]['mechs']['pas']['g'] = g[sec]
        netParams.cellParams['FoxP2_reduced']['secs'][sec]['mechs']['kdrin']['gkdrbar'] = gkdrbar[sec]
        netParams.cellParams['FoxP2_reduced']['secs'][sec]['mechs']['hin']['gbar'] = gbar[sec]
        netParams.cellParams['FoxP2_reduced']['secs'][sec]['mechs']['kapcb']['gkabar'] = gkabar[sec]
        netParams.cellParams['FoxP2_reduced']['secs'][sec]['mechs']['Nafx']['gnafbar'] = gnafbar[sec]
        netParams.cellParams['FoxP2_reduced']['secs'][sec]['mechs']['catcb']['gcatbar'] = gcatbar[sec]

    # For tuning the model and calibrate it
    for sec, secDict in netParams.cellParams['FoxP2_reduced']['secs'].items():
        if sec in cfg.tune:
            # vinit
            if 'vinit' in cfg.tune[sec]:
                secDict['vinit'] = cfg.tune[sec]['vinit']

            # mechs
            for mech in secDict['mechs']:
                if mech in cfg.tune[sec]:
                    for param in secDict['mechs'][mech]:
                        if param in cfg.tune[sec][mech]:
                            secDict['mechs'][mech][param] = cfg.tune[sec][mech][param]

            # geom
            for geomParam in secDict['geom']:
                if geomParam in cfg.tune[sec]:
                    secDict['geom'][geomParam] = cfg.tune[sec][geomParam]
else:
    cellRule = netParams.importCellParams(label='FoxP2_reduced', conds={'cellType':'FoxP2', 'cellModel':'HH_reduced'},
    fileName='cells/FS3.hoc', cellName='FoxP2', cellInstance=True)
    cellRule['secLists']['spiny'] = ['soma', 'dend']

#------------------------------------------------------------------------------
# Synaptic mechanism parameters
#------------------------------------------------------------------------------
netParams.synMechParams['NMDA'] = {'mod': 'MyExp2SynNMDABB', 'tau1NMDA': 2, 'tau2NMDA': 100, 'e': 0}
netParams.synMechParams['AMPA'] = {'mod': 'MyExp2SynBB', 'tau1': 0.5, 'tau2': 5*cfg.AMPATau2Factor, 'e': 0}

# ------------------------------------------------------------------------------
# VecStim inputs
# ------------------------------------------------------------------------------

if cfg.addVecStim:
    if 'Condition' in cfg.tune:
        cfg.Condition = cfg.tune['Condition']

    if 'RandomInitialization' in cfg.tune:
        cfg.Random = cfg.tune['RandomInitialization']

    ####
    # Load the spike trains
    if cfg.Go == 'Go':
        with open('cells/Results_Go.pkl', 'rb') as results:
            Results = pickle.load(results)
    elif cfg.Go == 'NoGo':
        with open('cells/Results_NoGo.pkl', 'rb') as results:
            Results = pickle.load(results)

    ###
    # Filter Results according to experimental data about connectivity
    minNumConn = cfg.RangeConnections[0]
    TrialList = [i for i in Results[cfg.Condition].keys()
                 if Results[cfg.Condition][i]['CellsPerTrial']>minNumConn]
    maxNumConn = cfg.RangeConnections[1]

    ###
    # Define population
    numTrials = len(TrialList)
    netParams.popParams['FoxP2'] = {'cellModel': 'HH_reduced', 'cellType': 'FoxP2', 'numCells': numTrials}

    for Trial in range(numTrials):
        cellsPerTrial = Results[cfg.Condition]['Trial_%d' % Trial]['CellsPerTrial']
        spkTimesIncre, cellIncre = [], []
        spkTimesDecre, cellDecre = [], []
        spkTimesNotChanging, cellNotChanging = [], []
        for keys, values in Results[cfg.Condition]['Trial_%d' % Trial].items():
            if keys != 'CellsPerTrial':
                cellID = Results[cfg.Condition]['Trial_%d' % Trial][keys]['ID']
                spikes = Results[cfg.Condition]['Trial_%d' % Trial][keys]['SimSpks']['Random_%d' % cfg.Random]
                if cellID == 'Increasing':
                    spkTimesIncre.append(spikes)
                    cellIncre.append(cellID)
                elif cellID == 'Decreasing':
                    spkTimesDecre.append(spikes)
                    cellDecre.append(cellID)
                elif cellID == 'Not Changing':
                    spkTimesNotChanging.append(spikes)
                    cellNotChanging.append(cellID)
        ####
        # Define the input populations
        netParams.popParams['Increasing_%d' % Trial] = {'cellModel': 'VecStim',
                                             'numCells': len(cellIncre),
                                             'spkTimes': spkTimesIncre}  # input from PT5B not changing firing rate
        netParams.popParams['Decreasing_%d' % Trial] = {'cellModel': 'VecStim',
                                             'numCells': len(cellDecre),
                                             'spkTimes': spkTimesDecre}  # input from PT5B not changing firing rate
        netParams.popParams['NotChanging_%d' % Trial] = {'cellModel': 'VecStim',
                                             'numCells': len(cellNotChanging),
                                             'spkTimes': spkTimesNotChanging}  # input from PT5B not changing firing rate

        Weights = {'Increasing' : cfg.AMPANMDAWeightsIncre, 'Decreasing' : cfg.AMPANMDAWeightsDecre,
                   'NotChanging' : cfg.AMPANMDAWeightsNotChanging}
    ####
    # Connect them
    PT5Bpops = [i for i in netParams.popParams.keys() if i is not 'FoxP2']

    for InputPops in PT5Bpops:
        connection = f'{InputPops}->FoxP2_%d' % int(InputPops.split('_')[1])
        numPreSyn = netParams.popParams[InputPops]['numCells']
        preSynSoma = np.sort(random.sample([i for i in range(numPreSyn)], int(cfg.somaProb*numPreSyn)))
        preSynDend = np.sort([i for i in range(numPreSyn) if i not in preSynSoma])

        netParams.connParams[connection+'_soma'] = {
                'preConds': {'popLabel': InputPops},
                'postConds': {'popLabel': 'FoxP2'},
                'weight': Weights[InputPops.split('_')[0]],
                'sec': 'soma',
                'connList': [[i,int(InputPops.split('_')[1])] for i in preSynSoma],
                'delay': cfg.delay,
                'loc': 0.5,
                'synMech': cfg.ESynMech}

        netParams.connParams[connection+'_dend'] = {
                'preConds': {'popLabel': InputPops},
                'postConds': {'popLabel': 'FoxP2'},
                'weight': Weights[InputPops.split('_')[0]],
                'sec': 'dend',
                'connList': [[i,int(InputPops.split('_')[1])] for i in preSynDend],
                'delay': cfg.delay,
                'loc': 'uniform(0,0.8)',
                'synMech': cfg.ESynMech}
else:
    netParams.popParams['FoxP2'] = {'cellModel': 'HH_reduced', 'cellType': 'FoxP2', 'numCells': 1}

#------------------------------------------------------------------------------
# Current inputs (IClamp)
#------------------------------------------------------------------------------
if cfg.addIClamp:
     for iclabel in [k for k in dir(cfg) if k.startswith('IClamp')]:
        ic = getattr(cfg, iclabel, None)  # get dict with params
        amps = ic['amp'] if isinstance(ic['amp'], list) else [ic['amp']]  # make amps a list if not already
        starts = ic['start'] if isinstance(ic['start'], list) else [ic['start']]  # make amps a list if not already
        for amp, start in zip(amps, starts):
            # add stim source
            netParams.stimSourceParams[iclabel+'_'+str(amp)] = {'type': 'IClamp', 'delay': start, 'dur': ic['dur'], 'amp': amp}
            # connect stim source to target
            netParams.stimTargetParams[iclabel+'_'+ic['pop']+'_'+str(amp)] = \
                {'source': iclabel+'_'+str(amp), 'conds': {'pop': ic['pop']}, 'sec': ic['sec'], 'loc': ic['loc']}
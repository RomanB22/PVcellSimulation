
"""
netParams.py

High-level specifications for M1 network model using NetPyNE
"""

from netpyne import specs
import pickle, json

netParams = specs.NetParams()   # object of class NetParams to store the network parameters

netParams.version = 49

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
netParams.sizeX = 100
netParams.sizeY = 1350  # cortical depth (will be converted to negative values)
netParams.sizeZ = 100


#------------------------------------------------------------------------------
# Cell parameters
#------------------------------------------------------------------------------
cellModels = ['HH_simple', 'HH_reduced', 'HH_full']
layer = {'2': [0.12,0.31], '4': [0.31,0.42], '5A': [0.42,0.52], '45A':[0.31,0.52], '5B': [0.52,0.77], '6': [0.77,1.0], 'long': [2.0,3.0]}  # normalized layer boundaries

#------------------------------------------------------------------------------
## Load cell rules previously saved using netpyne format
cellParamLabels = [] #['PV_reduced'] # list of cell rules to load from file
loadCellParams = cellParamLabels
saveCellParams = True #True

for ruleLabel in loadCellParams:
    netParams.loadCellParamsRule(label=ruleLabel, fileName='cells/'+ruleLabel+'_cellParams.pkl')

#------------------------------------------------------------------------------
# Specification of cell rules not previously loaded
# Includes importing from hoc template or python class, and setting additional params

#------------------------------------------------------------------------------
## PT5B full cell model params (700+ comps)
if 'PT5B_full' not in loadCellParams:
    ihMod2str = {'harnett': 1, 'kole': 2, 'migliore': 3}
    cellRule = netParams.importCellParams(label='PT5B_full', conds={'cellType': 'PT', 'cellModel': 'HH_full'},
      fileName='cells/PTcell.hoc', cellName='PTcell', cellArgs=[ihMod2str[cfg.ihModel], cfg.ihSlope], somaAtOrigin=True)
    nonSpiny = ['apic_0', 'apic_1']
    netParams.addCellParamsSecList(label='PT5B_full', secListName='perisom', somaDist=[0, 50])  # sections within 50 um of soma
    netParams.addCellParamsSecList(label='PT5B_full', secListName='below_soma', somaDistY=[-600, 0])  # sections within 0-300 um of soma
    for sec in nonSpiny: cellRule['secLists']['perisom'].remove(sec)
    cellRule['secLists']['alldend'] = [sec for sec in cellRule.secs if ('dend' in sec or 'apic' in sec)] # basal+apical
    cellRule['secLists']['apicdend'] = [sec for sec in cellRule.secs if ('apic' in sec)] # apical
    cellRule['secLists']['spiny'] = [sec for sec in cellRule['secLists']['alldend'] if sec not in nonSpiny]
    # Adapt ih params based on cfg param
    for secName in cellRule['secs']:
        for mechName,mech in cellRule['secs'][secName]['mechs'].items():
            if mechName in ['ih','h','h15', 'hd']:
                mech['gbar'] = [g*cfg.ihGbar for g in mech['gbar']] if isinstance(mech['gbar'],list) else mech['gbar']*cfg.ihGbar
                if cfg.ihModel == 'migliore':
                    mech['clk'] = cfg.ihlkc  # migliore's shunt current factor
                    mech['elk'] = cfg.ihlke  # migliore's shunt current reversal potential
                if secName.startswith('dend'):
                    mech['gbar'] *= cfg.ihGbarBasal  # modify ih conductance in soma+basal dendrites
                    mech['clk'] *= cfg.ihlkcBasal  # modify ih conductance in soma+basal dendrites
                if secName in cellRule['secLists']['below_soma']: #secName.startswith('dend'):
                    mech['clk'] *= cfg.ihlkcBelowSoma  # modify ih conductance in soma+basal dendrites
    # Reduce dend Na to avoid dend spikes (compensate properties by modifying axon params)
    for secName in cellRule['secLists']['alldend']:
        cellRule['secs'][secName]['mechs']['nax']['gbar'] = 0.0153130368342 * cfg.dendNa # 0.25
    cellRule['secs']['soma']['mechs']['nax']['gbar'] = 0.0153130368342  * cfg.somaNa
    cellRule['secs']['axon']['mechs']['nax']['gbar'] = 0.0153130368342  * cfg.axonNa # 11
    cellRule['secs']['axon']['geom']['Ra'] = 137.494564931 * cfg.axonRa # 0.005
    # Remove Na (TTX)
    if cfg.removeNa:
        for secName in cellRule['secs']: cellRule['secs'][secName]['mechs']['nax']['gbar'] = 0.0
    netParams.addCellParamsWeightNorm('PT5B_full', 'conn/PT5B_full_weightNorm.pkl', threshold=cfg.weightNormThreshold)  # load weight norm
    if saveCellParams:
        netParams.saveCellParamsRule(label='PT5B_full', fileName='cells/PT5B_full_cellParams.pkl')

#------------------------------------------------------------------------------
## PV cell params (3-comp)
if 'PV_reduced' not in loadCellParams:
    #cellRule = netParams.importCellParams(label='PV_reduced', conds={'cellType':'PV', 'cellModel':'HH_reduced'},
    #  fileName='cells/FS3.hoc', cellName='FScell1', cellInstance = True)
    cellRule = netParams.importCellParams(label='PV_reduced', conds={'cellType':'PV', 'cellModel':'HH_reduced'},
      fileName='cells/LTS3.hoc', cellName='FScell1', cellInstance = True)
    cellRule['secLists']['spiny'] = ['soma', 'dend']
    #netParams.addCellParamsWeightNorm('PV_reduced', 'conn/PV_reduced_weightNorm.pkl', threshold=cfg.weightNormThreshold)
    # cellRule['secs']['soma']['weightNorm'][0] *= 1.5
    #if saveCellParams: netParams.saveCellParamsRule(label='PV_reduced', fileName='cells/PV_reduced_cellParams.pkl')

for sec, secDict in netParams.cellParams['PV_reduced']['secs'].items():
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


Ra, cm, gkdrbar, gbar, gkabar, gnafbar, gcatbar, e, g = 226.8778027, \
1.4668784547221765, 0.074266696, 8.181467787372347e-07, 9.065582224086196e-06, \
0.6393160722589732, 0.00138129, -98.99999928, 4.65E-05

'''
cellRule= {conds: {cellType: 'PV', cellModel: 'HH_reduced'}, 
secs: {soma: {geom: {L: 42.0, nseg: 1, diam: 42.0, Ra: 226.8778027, cm: 1.4668784547221765}, 
topol: {}, mechs: {Nafx: {gnafbar: 0.6393160722589732, ar2: 1.0}, cadyn: {}, catcb: {gcatbar: 0.00138129}, 
hin: {K: 10.0, gbar: 8.181467787372347e-07, vhalf: -90.0}, kapcb: {gkabar: 9.065582224086196e-06}, 
kdrin: {gkdrbar: 0.074266696}, pas: {g: 4.65e-05, e: -98.99999928}}, 
ions: {ca: {e: 135.2142908191355, i: 5e-05, o: 2.0}, hi: {e: 0.0, i: 1.0, o: 1.0}, 
k: {e: -78.6041897101908, i: 54.4, o: 2.5}, na: {e: 50.0, i: 10.0, o: 140.0}}}, 
axon: {geom: {L: 113.22, nseg: 1, diam: 1.1, Ra: 150.0, cm: 1.2}, topol: {parentSec: 'soma', parentX: 0.5, childX: 0.0}, 
mechs: {Nafx: {gnafbar: 0.75, ar2: 1.0}, kdrin: {gkdrbar: 0.009}, pas: {g: 2.5e-05, e: -72.0750539439938}}, 
ions: {k: {e: -78.6041897101908, i: 54.4, o: 2.5}, na: {e: 50.0, i: 10.0, o: 140.0}}}, 
dend: {geom: {L: 176.0, nseg: 1, diam: 7.0, Ra: 150.0, cm: 1.2}, topol: {parentSec: 'soma', parentX: 0.0, childX: 0.0}, 
mechs: {Nafx: {gnafbar: 0.018, ar2: 1.0}, kapcb: {gkabar: 0.00875}, kdrin: {gkdrbar: 0.009}, 
pas: {g: 2.5e-05, e: -60.93996229584678}}, ions: {k: {e: -78.6041897101908, i: 54.4, o: 2.5}, 
na: {e: 50.0, i: 10.0, o: 140.0}}}}, secLists: {SectionList[0]: [], SectionList[1]: [], 
SectionList[16]: ['soma'], SectionList[17]: ['soma', 'axon', 'dend'], spiny: ['soma', 'dend']}, 
globals: {celsius: 23.0, hinf_catcb: 0.3068442636178239, ki0_k_ion: 140.0, ko0_k_ion: 3.82, minf_catcb: 0.0358812246203101}}
'''

#cellRule['secs']['soma']['vinit'] =  e
#cellRule['secs']['soma']['geom'] = {'cm': cm, 'Ra': Ra}
#cellRule['secs']['soma']['mechs']['pas'] = {'g': g, 'e': e}
#cellRule['secs']['soma']['mechs']['kdrin']['gkdrbar'] = gkdrbar
#cellRule['secs']['soma']['mechs']['hin']['gbar'] = gbar
#cellRule['secs']['soma']['mechs']['kapcb']['gkabar'] = gkabar
#cellRule['secs']['soma']['mechs']['Nafx']['gnafbar'] = gnafbar
#cellRule['secs']['soma']['mechs']['catcb']['gcatbar'] = gcatbar

#cellRule['secs']['axon']['vinit'] =  e
#cellRule['secs']['axon']['geom'] = {'cm': cm, 'Ra': Ra}
#cellRule['secs']['axon']['mechs']['pas'] = {'g': g, 'e': e}
#cellRule['secs']['dend']['vinit'] =  e
#cellRule['secs']['dend']['geom'] = {'cm': cm, 'Ra': Ra}
#cellRule['secs']['dend']['mechs']['pas'] = {'g': g, 'e': e}

#print(cellRule); #quit()

#------------------------------------------------------------------------------
# Population parameters
#------------------------------------------------------------------------------
#netParams.popParams['PT5B'] = {'cellModel': 'HH_full', 'cellType': 'PT', 'ynormRange': layer['5B'], 'numCells':1}
netParams.popParams['PV5B'] = {'cellModel': 'HH_reduced', 'cellType': 'PV', 'ynormRange': layer['5B'], 'numCells':1}

#------------------------------------------------------------------------------
# Synaptic mechanism parameters
#------------------------------------------------------------------------------
netParams.synMechParams['NMDA'] = {'mod': 'MyExp2SynNMDABB', 'tau1NMDA': 15, 'tau2NMDA': 150, 'e': 0}
netParams.synMechParams['AMPA'] = {'mod':'MyExp2SynBB', 'tau1': 0.05, 'tau2': 5.3*cfg.AMPATau2Factor, 'e': 0}
netParams.synMechParams['GABAB'] = {'mod':'MyExp2SynBB', 'tau1': 3.5, 'tau2': 260.9, 'e': -93}
netParams.synMechParams['GABAA'] = {'mod':'MyExp2SynBB', 'tau1': 0.07, 'tau2': 18.2, 'e': -80}
netParams.synMechParams['GABAASlow'] = {'mod': 'MyExp2SynBB','tau1': 2, 'tau2': 100, 'e': -80}

ESynMech = ['AMPA', 'NMDA']
SOMESynMech = ['GABAASlow','GABAB']
PVSynMech = ['GABAA']

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

#------------------------------------------------------------------------------
# NetStim inputs
#------------------------------------------------------------------------------
if cfg.addNetStim:
    for key in [k for k in dir(cfg) if k.startswith('NetStim')]:
        params = getattr(cfg, key, None)
        [pop, ynorm, sec, loc, synMech, synMechWeightFactor, start, interval, noise, number, weight, delay] = \
        [params[s] for s in ['pop', 'ynorm', 'sec', 'loc', 'synMech', 'synMechWeightFactor', 'start', 'interval', 'noise', 'number', 'weight', 'delay']]

        if synMech == ESynMech:
            wfrac = cfg.synWeightFractionEE
        elif synMech == SOMESynMech:
            wfrac = cfg.synWeightFractionSOME
        else:
            wfrac = [1.0]

        # add stim source
        netParams.stimSourceParams[key] = {'type': 'NetStim', 'start': start, 'interval': interval, 'noise': noise, 'number': number}

        # connect stim source to target
        # for i, syn in enumerate(synMech):
        netParams.stimTargetParams[key+'_'+pop] =  {
            'source': key,
            'conds': {'pop': pop, 'ynorm': ynorm},
            'sec': sec,
            'loc': loc,
            'synMech': synMech,
            'weight': weight,
            'synMechWeightFactor': synMechWeightFactor,
            'delay': delay}


#------------------------------------------------------------------------------
# VecStim inputs
#------------------------------------------------------------------------------

# #Set the trials
# # Create VecStims and connect them to the spinal cord interneuron
# trials_des=5
# #for index in range(len(next(iter(pop_decre_dict.values())))):  # Assuming all neurons have the same number of trials
# for index in range(trials_des):
#     for neuron_id in pop_decre_dict:
#         #print(pop_incre_dict[294][3])
#         #print(neuron_id)
#         trial_spikes = pop_decre_dict[neuron_id][index]
#         if trial_spikes.size>0:
#             vecstim_label = f'decre_neuron_{neuron_id}_trial_{index}'
#             netParams.stimSourceParams[vecstim_label] = {'type': 'VecStim', 'spikeTimes': trial_spikes}
#             netParams.stimTargetParams[f'{vecstim_label}->PV5B'] = {
#                 'source': vecstim_label,
#                 'conds': {'pop': 'PV5B'},
#                 'weight': 1.0,
#                 'delay': 1.0,
#                 'synMech': 'AMPA'
#             }
#     for neuron_id in pop_incre_dict:
#         print(neuron_id)
#         trial_spikes = pop_incre_dict[neuron_id][index]
#         if trial_spikes.size>0:
#             vecstim_label = f'incre_neuron_{neuron_id}_trial_{index}'
#             netParams.stimSourceParams[vecstim_label] = {'type': 'VecStim', 'spikeTimes': trial_spikes}
#             netParams.stimTargetParams[f'{vecstim_label}->PV5B'] = {
#                 'source': vecstim_label,
#                 'conds': {'pop': 'PV5B'},
#                 'weight': 1.0,
#                 'delay': 1.0,
#                 'synMech': 'AMPA'
#             }  s


with open('sim/InVivoFiringRate/CellPop.pkl', 'rb') as CellPop:
    neuron_properties = pickle.load(CellPop)

# Load the spike times
with open("sim/InVivoFiringRate/SpikeTimesMirrored.pkl","rb") as SpikeTimes:
    spike_times = pickle.load(SpikeTimes)

MinTrials=50

pop_decre = [i for i in neuron_properties.keys() if (neuron_properties[i]=='Increasing' and len(spike_times[i])>MinTrials)]
pop_incre = [i for i in neuron_properties.keys() if (neuron_properties[i]=='Decreasing' and len(spike_times[i])>MinTrials)]

pop_decre_dict = {}
pop_incre_dict = {}

#Populate the dictionaries based on labels
for key in spike_times:
    if key in pop_decre:
            pop_decre_dict[key] = spike_times[key]
    elif key in pop_incre:
            pop_incre_dict[key] = spike_times[key]

index=2
decre_spike_times_trials=[]
incre_spike_times_trials=[]


#for the desired index
for neuron in pop_decre_dict:
    neuron_lists=pop_decre_dict[neuron]
    decre_spike_times_trials.append(list(neuron_lists[index]))

for neuron in pop_incre_dict:
    neuron_lists=pop_incre_dict[neuron]
    incre_spike_times_trials.append(list(neuron_lists[index]))

decre_spike_times_trials = [[1000 + x for x in sublist] for sublist in decre_spike_times_trials]
incre_spike_times_trials = [[1000 + x for x in sublist] for sublist in incre_spike_times_trials]



if cfg.addVecStim:
    netParams.popParams['pop_incre'] = {'cellModel': 'VecStim', 'numCells': len(incre_spike_times_trials), 'spkTimes': incre_spike_times_trials}  # input from Ab_slow adapting type I

    netParams.connParams['pop_incre->PV5B'] = {
        'oneSynPerNetcon': True,
        'preConds': {'popLabel': 'pop_incre'}, 
        'postConds': {'popLabel': 'PV5B'},  
        'weight': 0.3221559,           
        'sec': 'soma',
        'probability': 0.9,
        'delay': 1.0,
        'loc': 0.5,
        'synMech': 'AMPA'} 
    
    netParams.popParams['pop_decre'] = {'cellModel': 'VecStim', 'numCells': len(decre_spike_times_trials), 'spkTimes': decre_spike_times_trials}  # input from Ab_slow adapting type I

    netParams.connParams['pop_decre->PV5B'] = {
        'oneSynPerNetcon': True,
        'preConds': {'popLabel': 'pop_decre'}, 
        'postConds': {'popLabel': 'PV5B'},  
        'weight': 0.3221559,           
        'sec': 'spiny',
        'probability': 0.9,
        'delay': 1.0,
        'loc': 0.5,
        'synMech': 'AMPA'} 


# netParams.connParams['pop_VecStim->PV5B'] = {
#     'preConds': {'popLabel': 'pop_VecStim'}, 
#     'postConds': {'popLabel': 'PV5B'},  
#     'weight': 0.5,           
#     'probability': 1.0,
#     'delay': 1.0,
#     'loc': 0.5,
#     'synMech': 'AMPA'} 


# spike_times_decre = [key for key in spike_times if key in pop_decre]
# spike_times_incre = [key for key in spike_times if key in pop_incre]

# netParams.popParams['pop_decre'] = {'cellModel': 'VecStim', 'numCells': len(pop_decre.keys()), 'spkTimes': cochlearSpkTimes}
# netParams.popParams['pop_decre'] = {'cellType': 'PYR', 'numCells': len(pop_decre)}
# netParams.popParams['pop_incre'] = {'cellType': 'PYR', 'numCells': len(pop_incre)}

# print(spkTimes_mirror[90][0])
# spkTimes_mirror[CellIndex][TrialIndex] = SpikeTimes list

# close files
SpikeTimes.close()
CellPop.close()


# 3 scenarios. Each one with all cells with at least 57 trials, and each scenario will have 57 trials.
# 1: Experimental: All cells with spike times. 57 Trials
# 2: Just the increasing population spike times (Experimental). "Ideal transmission". 57 Trials
# 3: Increasing Population (Experimental) + Mirrored decreasing population (Synthetic): "Worst case scenario". 57 Trials

# Caveat: Each cell has a different number of trials, so we need to look at a reasonable number of trials
# TODO: Select the neurons with at least 57 trials, and pick 57 trials as the number of trials.


#------------------------------------------------------------------------------
# Subcellular connectivity (synaptic distributions)
#------------------------------------------------------------------------------
if cfg.addSubConn:
    with open('conn/conn_dend_PT.json', 'rb') as fileObj: connDendPTData = json.load(fileObj)

    #------------------------------------------------------------------------------
    # Use subcellular distribution from L2/3 -> PT (Suter, 2015)
    lenY = 30
    spacing = 50
    gridY = list(range(0, -spacing*lenY, -spacing))
    synDens, _, fixedSomaY = connDendPTData['synDens'], connDendPTData['gridY'], connDendPTData['fixedSomaY']

    netParams.subConnParams['L2->PT'] = {
    'preConds': {'cellModel': 'NetStim'}, # all presyn inputs
    'postConds': {'pop': 'PT5B'},
    'sec': 'spiny',
    'groupSynMechs': ESynMech,
    'density': {'type': '1Dmap', 'gridX': None, 'gridY': gridY, 'gridValues': synDens['L2_PT'], 'fixedSomaY': fixedSomaY}}

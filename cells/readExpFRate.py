'''
    Code to read and save the experimental firing rates in nested dictionaries
    The format of the nested dictionary is:

    1. condition: InVivo_Go, InVivo_NoGo, OnlyIncre_Go, OnlyIncre_NoGo, MirrorDecre_Go, MirrorDecre_NoGo
        1.1. trial: Trial number according to experimental data. It changes from cell to cell
            1.1.1. cell: Cell number according to csv experimental files
                1.1.1.1. expFR: Experimental firing rate
                1.1.1.2. spktimes: Simulated spike times following an inhomogeneous Poisson process with probability given by expFR
                    1.1.1.2.1. Initialization: Random generation of spike trains, it could be more than one
                1.1.1.3 cell_id: 'Incre', 'Decre', 'Constant'. It depends on the variation of the average firing rate across trials
                                                              when comparing pre- and post-stimulus windows.
'''

import glob
import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from utils import *
import pickle as pkl


def LoadFiles(filename, preStim, postStim, sampFreq=31):
    cellNumber = [int(s) for s in re.findall(r'\d+', filename)]
    df = pd.read_csv(filename, sep=" ", index_col=None, header=0)
    # Find the CueTimes
    CueIndex = df.loc[df['cue'] == 1].index
    # Find the window around cue for the trial definition
    preWindowSize = int(preStim/sampFreq)
    postWindowSize = int(postStim/sampFreq)
    TrialWindowIdx = [[i-preWindowSize, i+postWindowSize] for i in CueIndex]
    time = [i*(1000/sampFreq) for i in range(preWindowSize+postWindowSize)]

    trialNumb = [[i]*len(time) for i in range(len(CueIndex))]
    trialNumb = [s for i in trialNumb for s in i]


    EpochDF = pd.concat([df.iloc[i[0]:i[1]] for i in TrialWindowIdx])
    EpochDF.insert(loc=0, column='Trial#', value=trialNumb)
    EpochDF.insert(loc=0, column='Time', value=time*len(CueIndex))
    EpochDF.insert(loc=0, column='Cell#', value=cellNumber*len(trialNumb))

    return EpochDF

def DefineCellID(df, preStim, postStim, plot=False):
    # It tells whether the neuron has 'Incre', 'Decre' or 'Constant' characteristics
    preStimRate = df.query('Time < %f' % preStim).groupby('Trial#').mean()['spikerate'].values
    postStimRate = df.query('Time >= %f' % postStim).groupby('Trial#').mean()['spikerate'].values

    rates = df.groupby('Time').mean()['spikerate'].values
    ratesStd = df.groupby('Time').std()['spikerate'].values

    minBtwPeriods = min([preStimRate.mean(), postStimRate.mean()])
    maxBtwPeriods = max([preStimRate.mean(), postStimRate.mean()])
    semiDiffBtwPeriods = (postStimRate.mean() - preStimRate.mean())/2

    if maxBtwPeriods/minBtwPeriods < 1.1: # If max value is at least 10% higher than min
        df.insert(loc=0, column='Cell_ID', value=['Not Changing'] * len(df.index))
    else:
        if semiDiffBtwPeriods > 0:
            df.insert(loc=0, column='Cell_ID', value=['Increasing'] * len(df.index))
        elif semiDiffBtwPeriods < 0:
            df.insert(loc=0, column='Cell_ID', value=['Decreasing'] * len(df.index))

    if plot:
        for i in np.unique(df['Trial#']):
            plt.plot(df.loc[df['Trial#'] == i]['Time'], df.loc[df['Trial#'] == i]['spikerate'], color='gray', alpha=0.5)

        yerr0 = rates-ratesStd
        yerr0[yerr0<0] = 0
        plt.fill_between(df.loc[df['Trial#'] == 0]['Time'], yerr0, rates+ratesStd, color='tab:blue', alpha=0.5)
        plt.plot(df.loc[df['Trial#'] == 0]['Time'], rates, color='tab:blue', label = 'Average', linewidth = 5)
        plt.hlines(preStimRate.mean(),0,preStim, color='r')
        plt.hlines(postStimRate.mean(), preStim, preStim+postStim, color='r')
        plt.legend()

        plt.savefig("Average_Cell%d.png" % np.unique(df["Cell#"])[0])
        plt.clf()

    return df

def SimulateSpikes(ExpFR, MirroredExpFR, time, CellPop, condition):
    if (condition == 'MirrorDecre_Go' or condition == 'MirrorDecre_NoGo') and CellPop in ['Decreasing']:
        rates = MirroredExpFR
    elif (condition == 'OnlyIncre_Go' or condition == 'OnlyIncre_NoGo') and CellPop not in ['Increasing']:
        rates = np.zeros(np.shape(ExpFR))
    else:
        rates = ExpFR

    tstop = time[-1]
    spkTimes = [x for x in inh_poisson_generator(rates, time, tstop)]

    return spkTimes


def CreateDict(df, preStim, numInitializations=10, Conditions = ('InVivo_Go', 'OnlyIncre_Go', 'MirrorDecre_Go')):
    # simulate spike trains
    # df has columns ['Cell_ID', 'Cell#', 'Time', 'Trial#', 'spikerate', 'ITI', 'cue', 'push',
    #        'reward']
    dictSpikes = {}
    dictSpikes['Time'] = np.unique(df['Time'].values)
    for cond in Conditions:
        dictSpikes[cond] = {}
        for trial in np.unique(df['Trial#']):
            cellsId = np.unique(df.loc[df['Trial#'] == trial]['Cell#'])
            numCellPerTrial = len(cellsId)
            dictSpikes[cond]['Trial_%d' % int(trial)] = {}
            dictSpikes[cond]['Trial_%d' % int(trial)]['CellsPerTrial'] = numCellPerTrial
            for cell in cellsId:
                dictSpikes[cond]['Trial_%d' % int(trial)]['Cell_%d' % int(cell)] = {}
                CellPop = np.unique(df.loc[df['Cell#'] == int(cell)]['Cell_ID'])[0]
                dictSpikes[cond]['Trial_%d' % int(trial)]['Cell_%d' % int(cell)]['ID'] = CellPop
                # Save ExpFR and mirrored firing rate
                ExpFR = df.loc[(df['Cell#'] == int(cell)) & (df['Trial#'] == trial)]['spikerate'].values
                PreStimFR = \
                    df.query('Time <= %f' % preStim).loc[(df['Cell#'] == int(cell)) & (df['Trial#'] == trial)]['spikerate'].values
                MirroredExpFR = np.concatenate((PreStimFR, np.flip(PreStimFR)), axis=0)
                if len(MirroredExpFR) < len(ExpFR):
                    MirroredExpFR = np.pad(MirroredExpFR,(0,len(ExpFR) - len(MirroredExpFR)), 'constant')
                elif len(MirroredExpFR) > len(ExpFR):
                    MirroredExpFR = MirroredExpFR[:len(ExpFR)]
                ###
                dictSpikes[cond]['Trial_%d' % int(trial)]['Cell_%d' % int(cell)]['ExpFR'] = ExpFR
                dictSpikes[cond]['Trial_%d' % int(trial)]['Cell_%d' % int(cell)]['MirroredExpFR'] = MirroredExpFR
                dictSpikes[cond]['Trial_%d' % int(trial)]['Cell_%d' % int(cell)]['SimSpks'] = {}
                for i in range(numInitializations): # We can parallelize this step in the future. No need to do it in a for loop
                    dictSpikes[cond]['Trial_%d' % int(trial)]['Cell_%d' % int(cell)]['SimSpks']['Random_%d' % i] = \
                    SimulateSpikes(ExpFR, MirroredExpFR, np.unique(df['Time'].values), CellPop, cond)

    return dictSpikes

if __name__=="__main__":
    folder = 'Go'
    Conditions = ('InVivo_%s' % folder, 'OnlyIncre_%s' % folder, 'MirrorDecre_%s' % folder)
    path = r'CSN_spike_data_%s/' % folder  # use your path
    all_files = glob.glob(os.path.join(path, "*.txt"))
    preStim, postStim = 1800, 1800

    df = pd.DataFrame()
    for filename in all_files:
        dfAux = LoadFiles(filename, preStim, postStim)
        dfAux = DefineCellID(dfAux, preStim, postStim, plot=False)
        df = pd.concat([df,dfAux], ignore_index=True)

    Results = CreateDict(df, preStim, numInitializations=1, Conditions=Conditions)
    print(Results[Conditions[0]].keys())
    print(Results[Conditions[0]]['Trial_1'].keys())
    print(Results[Conditions[0]]['Trial_1']['CellsPerTrial'])

    with open('Results_%s.pkl' % folder, 'wb') as f:
        pkl.dump(Results, f)




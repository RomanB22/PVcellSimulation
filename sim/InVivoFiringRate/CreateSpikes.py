from utils import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle as pkl

df = pd.read_csv('RatesAllCells.csv')

spkTimes_df = {}
preStim=1000
postStim=1000
pushTime = []

CellPop = {}
for cell in np.unique(df['Cell#']):
    cell = int(cell)
    cellDF = df.loc[df['Cell#']==cell]
    rates = cellDF['spikerate'].values

    meanFR = []

    for trial in np.unique(cellDF['Trial#']):
        trial = int(trial)
        cellDF2 = cellDF.loc[cellDF['Trial#'] == trial]
        rates2 = cellDF2['spikerate'].values
        times2 = cellDF2['Time'].values
        timeMin = np.min(times2)
        cellDF3 = cellDF2.loc[cellDF2['cue'] == 1]
        CueTime = cellDF3['Time'].values
        index = np.argwhere((times2 >= CueTime - preStim) & (times2 < CueTime+ postStim)).flatten()
        meanFR.append(rates2[index])

    minIndex = np.min([len(i) for i in meanFR])
    meanFR2 = [i[:minIndex] for i in meanFR]
    meanFR2 = np.mean(meanFR2,axis=0)
    stdFR2 = np.std([i[:minIndex] for i in meanFR],axis=0, ddof=1)/np.sqrt(trial+1)
    yerr0 = meanFR2 - stdFR2
    yerr0[yerr0<0] = 0
    yerr1 = meanFR2 + stdFR2

    timeBeforeCue = (times2[index[:minIndex]]-CueTime)<=0
    timeAfterCue = (times2[index[:minIndex]]-CueTime)>0
    meanBeforeCue = np.mean(meanFR2[timeBeforeCue])
    meanAfterCue = np.mean(meanFR2[timeAfterCue])

    semiDiffBtwPeriods = (meanAfterCue - meanBeforeCue)/2
    minBtwPeriods = min([meanAfterCue, meanBeforeCue])
    maxBtwPeriods = max([meanAfterCue, meanBeforeCue])
    meanBtwPeriods = (meanAfterCue + meanBeforeCue)/2

    if abs(meanBtwPeriods)/maxBtwPeriods >= 0.9: # If mean and max value are similar in at least 90%
        plt.title("Not correlated")
        CellPop[cell] = "Not Correlated"
    else:
        if semiDiffBtwPeriods > 0:
            plt.title("Increasing")
            CellPop[cell] = "Increasing"
        elif semiDiffBtwPeriods < 0:
            plt.title("Decreasing")
            CellPop[cell] = "Decreasing"

    plottingTimes= times2[index[:minIndex]]-CueTime

    plt.axvline(x=0,ls='--',c='k')
    plt.fill_between(plottingTimes, yerr0, yerr1, color='tab:blue', alpha=0.5)
    plt.plot(plottingTimes, meanFR2, color='tab:blue', label=r'$mean \pm std_m$')
    plt.plot(plottingTimes[timeBeforeCue], meanBeforeCue*np.ones(np.shape(plottingTimes[timeBeforeCue])), color='r',label='period mean value')
    plt.plot(plottingTimes[timeAfterCue], meanAfterCue*np.ones(np.shape(plottingTimes[timeAfterCue])), color='r')
    plt.xlabel("Time (ms)")
    plt.ylabel("Firing Rate (Hz)")
    plt.legend()
    ax = plt.gca()
    ax.spines[['right', 'top']].set_visible(False)
    plt.savefig("PlotsCells/AverageFiringRate_cell_%d.png" % cell)
    plt.close()

# create a binary pickle file
f = open("CellPop.pkl","wb")

# write the python object (dict) to pickle file
pkl.dump(CellPop,f)

# close file
f.close()

for cell in np.unique(df['Cell#']):
    cell = int(cell)
    spkTimes_df[cell] = {}
    cellDF = df.loc[df['Cell#']==cell]
    rates = cellDF['spikerate'].values
    times = cellDF['Time'].values
    tstop = times[-1]
    spkTimes = [x for x in inh_poisson_generator(rates, times, tstop)]

    #plt.plot(times, rates,'--k')
    #plt.scatter(times[cellDF['push']==1], rates[cellDF['push']==1])
    #plt.scatter(spkTimes, np.ones(len(spkTimes)))
    #plt.ylim([0, max(rates)])
    cellDF4 = cellDF.loc[cellDF['cue'] == 1]
    #CueTime = cellDF4['Time'].values
    #for xc in CueTime:
    #    plt.axvline(x=xc)
    #plt.savefig('PlotsCells/Cell_%d_FullTrial.png' % cell)
    #plt.close()
    for trial in np.unique(cellDF['Trial#']):
        trial = int(trial)
        cellDF2 = cellDF.loc[cellDF['Trial#'] == trial]
        rates2 = cellDF2['spikerate'].values
        times2 = cellDF2['Time'].values
        tstop2 = times2[-1]

        cellDF3 = cellDF2.loc[cellDF2['cue'] == 1]
        CueTime = cellDF3['Time'].values
        pushTime.append(max(times2[cellDF2['push'] == 1])-CueTime)

        # Two ways of generating spikes. For whole session and then selecting the ones of interest
        # Vs. creating the poisson for each trial (gives more spikes the second way)
        Method=2
        if Method == 2:
            spkTimes2 = [x for x in inh_poisson_generator(rates2, times2, tstop2)]
            index = np.argwhere((spkTimes2 >= CueTime - preStim) & (spkTimes2 < CueTime+ postStim)).flatten()
            spkTimes_df[cell][trial] = [float(spkTimes2[i]) for i in index]-CueTime
        else:
            #index = np.argwhere((spkTimes >= min(times2)) & (spkTimes < max(times2))).flatten()
            index = np.argwhere((spkTimes >= CueTime-preStim) & (spkTimes < CueTime+postStim) ).flatten()

            if len(index) > 0:
                minidx = int(min(index))
                maxidx = int(max(index))
                spkTimes_df[cell][trial] = spkTimes[minidx:maxidx ]-CueTime
                spkTimes2 = spkTimes[minidx:maxidx]
            else:
                spkTimes2 = []
        # plt.plot(times2-min(times2), rates2,'--k')
        # plt.scatter(times2[cellDF2['push'] == 1]-min(times2), rates2[cellDF2['push'] == 1])
        # plt.scatter(spkTimes2-min(times2),np.ones(len(spkTimes2)))
        # for xc in CueTime:
        #     plt.axvline(x=xc-min(times2))
        # plt.ylim([0, max(rates)])
        # plt.xlim([CueTime-1000-min(times2), CueTime-min(times2)+1000])
        # plt.savefig('PlotsCells/Cell_%d_Trial_%d.png' % (cell,trial) )
        # plt.close()

for cell in spkTimes_df.keys():
    keys = list(spkTimes_df[cell].keys())
    for key in keys:
        print(cell, key, spkTimes_df[cell][key])

print(spkTimes_df)

# create a binary pickle file
f = open("SpikeTimes.pkl","wb")

# write the python object (dict) to pickle file
pkl.dump(spkTimes_df,f)

# close file
f.close()
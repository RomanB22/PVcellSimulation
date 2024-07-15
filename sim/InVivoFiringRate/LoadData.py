import pandas as pd
import glob
import os
import re
import numpy as np

path = r'CSN_spike_data/' # use your path
all_files = glob.glob(os.path.join(path, "*.txt"))

li = []
minIndex = []
for filename in all_files:
    cellNumber=[int(s) for s in re.findall(r'\d+', filename)]
    df = pd.read_csv(filename, sep=" ", index_col=None, header=0)
    minIndex.append(len(df.index))
    # Find the trial number
    EndTrial = df.loc[(df['reward'].diff()==-1) & (df['ITI'].diff()==1)].index
    EndTrialAux = np.append(np.array([0]),EndTrial)
    EndTrialAux = np.append(EndTrialAux,np.array(df.index[-1]+1))
    trialNumb = np.zeros(len(df.index),dtype=int)
    trialIdx=0
    sampFreq=31 # 31 Hz
    timeWindow=1000/sampFreq
    time = np.array([i*timeWindow for i in range(len(df.index))])
    for i in range(len(EndTrialAux)-1):
        trialNumb[EndTrialAux[i]:EndTrialAux[i+1]]=trialIdx
        trialIdx+=1

    df.insert(loc=0, column='Trial#', value=trialNumb)
    df.insert(loc=0, column='Time', value=time)
    df.insert(loc=0, column='Cell#', value=cellNumber*np.ones(len(df.index),dtype=int))
    li.append(df)
print(min(minIndex), max(minIndex))
frame = pd.concat(li, axis=0, ignore_index=True)

frame.to_csv('RatesAllCells.csv', index=False)

import midas.file_reader
from datetime import datetime
import numpy as np
import cygno as cy
from tqdm.notebook import tqdm
import pandas as pd
from pmt_lib import *
from scipy.signal import find_peaks


def moving_average(array, window_size = 5, drop_NaN = True):
    tmp = pd.Series(array)
    moving_average = (tmp.rolling(window=window_size).mean().to_numpy())
    
    if drop_NaN:
        moving_average = (pd.Series(moving_average).dropna()).to_numpy()
    
    return moving_average


def moving_average_4ch(arrays, window_size = 5, channels = 4, drop_NaN = True):
    
    ma = []
    for i in range(channels):
        ma.append(moving_average(arrays[i], window_size, drop_NaN=drop_NaN))
        
    return np.array(ma)


def four_channel_peakfinder(wfs, width = 5, prominence = True):
    
    
    if prominence:
        peaks1, peaks_dict1 = find_peaks(wfs[0], width=width, prominence=max(wfs[0])/4)
        peaks2, peaks_dict2 = find_peaks(wfs[1], width=width, prominence=max(wfs[1])/4)
        peaks3, peaks_dict3 = find_peaks(wfs[2], width=width, prominence=max(wfs[2])/4)
        peaks4, peaks_dict4 = find_peaks(wfs[3], width=width, prominence=max(wfs[3])/4)
    else: 
        peaks1, peaks_dict1 = find_peaks(wfs[0], width=width, prominence=0.01)
        peaks2, peaks_dict2 = find_peaks(wfs[1], width=width, prominence=0.01)
        peaks3, peaks_dict3 = find_peaks(wfs[2], width=width, prominence=0.01)
        peaks4, peaks_dict4 = find_peaks(wfs[3], width=width, prominence=0.01)

    
    return [peaks1, peaks2, peaks3, peaks4], [peaks_dict1, peaks_dict2, peaks_dict3, peaks_dict4]

def find_majority2_peaks(pk, window=15):
    multi=[] 

    for index in range(len(pk)): #loop over the 4 arrays

        pk_tmp = pk[:index]+ pk[index+1:] #create the list of the other 3 arrays

        for p in pk[index]:               #loop over the elements of the selectes array
            cont=-1
            for pp in pk_tmp:             #loop over the arrays not selected
                for p2 in pp:             #loop over the value of the array
                    flag=False
                    if cont==1: break     #break if the majority2 was already found
                    if np.abs(p-p2)<15:   #condition: if the peaks of 2 different channels
                        cont=1              # are in a 15 sample window there is a match
                        
                        if len(multi)!=0:
                            not_found=-1
                            for m in multi:          #need to check if that peak was already evaluated
                                if np.abs(p-m)<15:      #if there is a value in the window break
                                    flag=True
                                    break
                                if flag==True:
                                    break

                            if flag==False: multi.append(p)  #if no value in the window was found, append it

                        else: multi.append(p)
    multi = np.sort(np.array(multi))
    return multi


def run_peak_analizer(WFs, run=0, av_size=5, prominence=True, window=15):
    r = 1/4096
    
    pictures = np.arange(len(WFs))
    majority2_peaks = []
    run_arr=np.array([], dtype=int)
    pic_arr=np.array([], dtype=int)
    trg_arr=np.array([], dtype=int)

    for pic in tqdm(pictures): # LOOP OVER THE PICTURES
        triggers = np.arange(len(WFs[pic]))

        for trg in triggers: # LOOP OVER THE TRIGGERS
            # find the peaks
            wfs = WFs[pic][trg][1:5]
            baseline = np.array([np.mean(wfs[0][:100]), np.mean(wfs[1][:100]), 
                                 np.mean(wfs[2][:100]), np.mean(wfs[3][:100])])
            rms = np.array([np.std(wfs[0][:100]), np.std(wfs[1][:100]),
                            np.std(wfs[2][:100]), np.std(wfs[3][:100])])
            for i in range(4):
                wfs[i] = wfs[i]-baseline[i]

            wf_averaged = moving_average_4ch(wfs, window_size=av_size)
            peaks, peaks_dict = four_channel_peakfinder(-r*wf_averaged, width=5, prominence=prominence)

            # CHECK MAJORITY 2 peaks
            majority2_peaks.append(find_majority2_peaks(peaks, window=window))
            
            run_arr=np.append(run_arr, run)
            pic_arr=np.append(pic_arr, pic)
            trg_arr=np.append(trg_arr, trg)

    df={'run':run_arr, 'picture':pic_arr.astype(int), 'trigger':trg_arr.astype(int), 'peaks':peaks, 'mj2 peaks':majority2_peaks}
    df_int=pd.DataFrame(df)
    return df_int
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 28 15:27:15 2021

@author: Luis González
"""
import pandas as pd
import tkinter as tk
from tkinter import filedialog
import numpy as np
from scipy.fft import fft, ifft
import matplotlib.pyplot as plt
import cmath
import math
import os
from functions import RotatingCoilAnalysisTurn



print('Pick Raw Data file')
# file = filedialog.askopenfilename()
file=r'C:\Users\Luis González\cernbox\Work\CIEMAT-\MCBXF\FFMM\Python Post-processing\Matlab\data for luis\MCBXFB-P2-inner-5A-ctr-pos15_Run_01_I_-5.00A\MCBXFB-P2-inner-5A-ctr-pos15_Run_01_I_-5.00A_raw_measurement_data.txt'
file_name=file.split('/')[-1]
file_folder = '\\'.join(file.split('\\')[0:-1])

print('Pick Sensitivities file')
# sens = filedialog.askopenfilename()
sens=r'C:\Users\Luis González\cernbox\Work\CIEMAT-\MCBXF\FFMM\Python Post-processing\Matlab\data for luis\Kn_R45_130_N1_0001_A_AC.txt'
sens_name=file.split('\\')[0:-1]
sens_folder = '\\'.join(file.split('\\')[0:-1])
p_turn=512
num_FDI=2
MagOrder=1
Rref=0.050
AnalysisOptions = 'dri rot nor cel'

sensitivities = pd.read_csv(sens, sep='\s+',names=['abs_real','abs_im','cmp_real','cmp_im'])
knAbs=np.array(sensitivities['abs_real']+ 1j * sensitivities['abs_im'])
knCmp=np.array(sensitivities['cmp_real']+ 1j * sensitivities['cmp_im'])

#def ContinuousRotatingCoilAnalysis (fileName, OutputFolder, p_turn, knAbs, knCmp, num_FDI , MagOrder, Rref, AnalysisOptions):

df = pd.read_csv(file, sep='\s+',names=['Time','df_abs','df_cmp','curr'])
turns=len(df)/p_turn
turnlst=[]
i=0

for turn in np.arange(int(turns)):

    turn_data=df.iloc[i*p_turn:(i+1)*p_turn].reset_index()
    turnlst.append(turn_data)
    i=i+1
# take the turn that we want, e.g. 2
# turn=turndict[9]





NAbs=pd.DataFrame()
SAbs=pd.DataFrame()
NCmp=pd.DataFrame()
SCmp=pd.DataFrame()
Position=pd.DataFrame(columns=("Turn","X","Y","Phi"))

#Analise the signal for each turn


i=0 
for turn in turnlst:
    #drift correction----------------------------------------------------------------------------
    if 'dri'in AnalysisOptions:
        Fabs = (turn['df_abs']-turn['df_abs'].mean()).cumsum()-turn['df_abs'].cumsum().mean()
        Fcmp = (turn['df_cmp']-turn['df_cmp'].mean()).cumsum()-turn['df_cmp'].cumsum().mean()
    else:
        Fabs=turn['df_abs'].cumsum()
        Fcmp=turn['df_cmp'].cumsum()
    #Call the analysis funtion of rotatin coil----------------------------------------------------------------------------
    [Norm_abs,
     Skew_abs,
     Norm_cmp,
     Skew_cmp,
     x,
     y,
     PhiOut]=RotatingCoilAnalysisTurn(Fabs, Fcmp, knAbs, knCmp, MagOrder, Rref, AnalysisOptions)
    
    # Add the Multipoles corresponding to each turn to different Dataframes depending on their NOrmal/Skew and Absolute/Compensated
    
    NAbs["Turn_"+str(i)]=Norm_abs.tolist()
    SAbs["Turn_"+str(i)]=Skew_abs.tolist()
    NCmp["Turn_"+str(i)]=Norm_cmp.tolist()
    SCmp["Turn_"+str(i)]=Skew_cmp.tolist()
    
    Position.loc[len(Position)]=[i+1,x,y,PhiOut]
    
    i=i+1

Multip_list=[NAbs,SAbs,NCmp,SCmp] # Add number of the Multipole
for dfr in Multip_list:
    
    dfr=dfr.insert(0,"n",dfr.index+1)

Multip_list=[NAbs,SAbs,NCmp,SCmp] # Plot multiples Obtained at each turn
for dfr in Multip_list:
    plt.figure()
    for column in dfr.columns[1:]:
        plt.plot(dfr["n"],dfr[column])
        
# Save the dataframes
NAbs.to_csv(file.replace("raw_measurement_data","C_Normal_Absolute"))
SAbs.to_csv(file.replace("raw_measurement_data","C_Skew_Absolute"))
NCmp.to_csv(file.replace("raw_measurement_data","C_Normal_Compensated"))
SCmp.to_csv(file.replace("raw_measurement_data","C_Skew_Compensated"))
Position.to_csv(file.replace("raw_measurement_data","Position"))



    






    
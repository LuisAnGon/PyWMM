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
from functions import ReadRoxie, RotatingCoilAnalysisTurn,Parameters, ContinuousRotatingCoilAnalysis, GetSensitivities,ReadFromC,Differences



step=65

#folder=r'C:\Users\Luis González\cernbox\Work\CIEMAT-\MCBXF\FFMM\rotating coil processing\Python\Data at 927\Scan-Fringe'
folder=r'C:\Single Coils Centre Rot90\ID_LC_90_ctr'


senspath=folder+'\\Kn_R45_130_N1_0001_A_AC.txt'    
[knAbs,knCmp]=GetSensitivities(senspath)

ev_mul=pd.DataFrame()

for Phi in np.linspace(-0.2028+(0*np.pi),-0.2028+(0*np.pi),1):
# for Phi in np.array([1.3721+(0*np.pi),
#                       1.3721+(0.5*np.pi),
#                       1.3721+(1*np.pi),
#                       1.3721+(1.5*np.pi)]):
# for Phi in np.linspace((0)*np.pi,(2)*np.pi,200)

    meas=[] #Each folder Corresponds to one measurement in one position
    for name in os.listdir(folder):
        if os.path.isdir(os.path.join(folder, name)):
            meas.append(name)
    
    
    Av_pos=pd.DataFrame()
    Av_neg=pd.DataFrame()
    
    
    i_pos=1
    paso=[]
    for measurement in meas:  # DIFFERENT MEASUREMENTS Meanining measurements at different positions
        paso.append(step*(i_pos-1))
        
        for item in os.listdir(folder+'\\'+measurement): #Get the parameters from the exported file from C++
            if item.split('_')[-1]=='Parameters.txt':# Read the parameters at each position
                parameterpath=folder+'\\'+measurement+'\\'+item
                [p_turn, num_FDI,MagOrder,Rref,AnalysisOptions]=Parameters(parameterpath)
        
        
        for dire in [x[0] for x in os.walk(folder+'\\'+measurement)][1:]:# We go to each position
            
            Current=dire.split('_')[-1] # Read the current at which the measurement was done
           
            for doc in os.listdir(dire): #We obtain every measurment done at each position
                
                if 'raw_measurement_data' in doc:
                    datafile=dire+'\\'+doc         
                    [PhiOut,Av]=ContinuousRotatingCoilAnalysis (datafile, p_turn, knAbs, knCmp, MagOrder, Rref, AnalysisOptions,step,Phi) #Get the average of all the Rotation coil turns
                             
        
            Av=Av.rename(columns={"Average":i_pos-1}) #Rename each column with the number of position
            
            #----------------------------- Put it in its correct dataframe according to the current of the magnet
            
            if float(Current.split('A')[0])>0:
                Av_pos=pd.concat([Av_pos,Av],axis=1)
                # Av_skew_pos=Av_skew_pos.append(SCmp_ave)
                
            if float(Current.split('A')[0])<0: 
                Av_neg=pd.concat([Av_neg,Av],axis=1)
                # Av_skew_neg=Av_skew_neg.append(SCmp_ave)    
        i_pos=i_pos+1     
    
    # Specify the position of each measurement according to the number of position
    Av_neg.loc["position"]=paso
    Av_pos.loc["position"]=paso
    
    Av=pd.DataFrame()
    for col in Av_neg.columns:
        Av[col]=(Av_pos[col]+Av_neg[col])/2
        
    
    # Create DataFrames with non normalized multipoles
    Av_NN=pd.DataFrame()
    Av_NS=pd.DataFrame()
    
    
    # print(Av)
    
    for row in Av.iterrows():
        if row[0]=='position':
            Av_NN['position']=Av.loc[row[0]]
            Av_NS['position']=Av.loc[row[0]]
            
        if row[0][0]=='B':
            Av_NN['B1']=Av.loc[row[0]]        
        elif row[0][0]=='b':
            Av_NN[row[0]]=Av.loc[row[0]]        
        elif row[0][0]=='A':
            Av_NS['A1']=Av.loc[row[0]]        
        elif row[0][0]=='a':
            Av_NS[row[0]]=Av.loc[row[0]]
    Av_NN.to_excel(folder+'\\Av_NN'+str(Phi.round(4))+'pi.xlsx')
    Av_NS.to_excel(folder+'\\Av_NS'+str(Phi.round(4))+'pi.xlsx')
    
    
    [RoxieMP,Roxieall]=ReadRoxie(folder)
    
    # print('ROXIE ALL:',RoxieMP)
    if sum(RoxieMP.index.isin(['A16']))>0:
        RoxieMP=RoxieMP.drop('A16',axis=0)
    RoxieMP.index=["B1","b2","b3","b4","b5","b6","b7","b8","b9","b10","b11","b12","b13","b14","b15",
                   "A1","a2","a3","a4","a5","a6","a7","a8","a9","a10","a11","a12","a13","a14","a15"]
    pos=Av.loc[["B1","b2","b3","b4","b5","b6","b7","b8","b9","b10","b11","b12","b13","b14","b15",
                "A1","a2","a3","a4","a5","a6","a7","a8","a9","a10","a11","a12","a13","a14","a15"]][Av.columns[int(len(Av.columns)/2)]]
       
        
    comp_Roxie_Meas=pd.concat([RoxieMP,pos],axis=1)
    comp_Roxie_Meas=comp_Roxie_Meas.round(6)
    comp_Roxie_Meas.columns=['Roxie','MM']
        
        
        
    comp_Roxie_Meas_B=comp_Roxie_Meas.loc[["b2","b3","b4","b5","b6","b7","b8","b9","b10","b11","b12","b13","b14","b15"]]
    comp_Roxie_Meas_B['n']=np.arange(2,16)
    
    comp_Roxie_Meas_A=comp_Roxie_Meas.loc[["a2","a3","a4","a5","a6","a7","a8","a9","a10","a11","a12","a13","a14","a15"]]
    comp_Roxie_Meas_A['n']=np.arange(2,16)
    
    plt.figure()
    plt.bar(comp_Roxie_Meas_B['n'],comp_Roxie_Meas_B['MM'], alpha=0.5,label='MM')
    # plt.bar(comp_Roxie_Meas_B['n'],comp_Roxie_Meas_B['Roxie'], alpha=0.5,label='Roxie')
    plt.title('Normal Multipoles: Comparison MM and Roxie - Rot Angle='+str((PhiOut/np.pi).round(5))+'pi')
    plt.xlabel('n')
    plt.ylabel('b (units)')
    plt.legend()
    
    plt.figure()
    plt.bar(comp_Roxie_Meas_A['n'],comp_Roxie_Meas_A['MM'], alpha=0.5,label='MM')
    # plt.bar(comp_Roxie_Meas_A['n'],comp_Roxie_Meas_A['Roxie'], alpha=0.5,label='Roxie')
    plt.title('Skew Multipoles: Comparison MM and Roxie - Rot Angle='+str((PhiOut/np.pi).round(5))+'pi')
    plt.xlabel('n')
    plt.ylabel('a (units)')
    plt.legend()
    












                



    






    
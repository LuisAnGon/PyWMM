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
folder=r'C:\Single Coils Rot 0\ID_UC_Nom_ctr'


senspath=folder+'\\Kn_R45_130_N1_0001_A_AC.txt'    
[knAbs,knCmp]=GetSensitivities(senspath)

ev_mul=pd.DataFrame()


a=0
mul=0
Angle=[]
Diff1=[]
Diff2=[]
Diff3=[]
Diff4=[]
Diff5=[]
Diff6=[]
Diff7=[]
Diff8=[]
Diff9=[]
Diff10=[]
Diff11=[]
Diff12=[]
Diff13=[]
Diff14=[]
Diff15=[]
Diff16=[]
while a==0:
    Phi=mul*np.pi
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
    
    # Av_NN.to_excel(folder+'\\Av_NN'+str(Phi)+'pi.xlsx')
    # Av_NS.to_excel(folder+'\\Av_NS'+str(Phi)+'pi.xlsx')
     
    # dif1=(Av_NN['b3'].iloc[0].round(5))-(Av_NS['a3'].iloc[0].round(5))

    
    Angle.append(mul)
    # Diff1.append(dif1)
    Diff1.append(Av_NN['b2'].iloc[0].round(5))
    Diff2.append(Av_NN['b3'].iloc[0].round(5))
    Diff3.append(Av_NN['b4'].iloc[0].round(5))
    Diff4.append(Av_NN['b5'].iloc[0].round(5))
    Diff5.append(Av_NN['b6'].iloc[0].round(5))
    Diff6.append(Av_NN['b7'].iloc[0].round(5))
    Diff7.append(Av_NN['b8'].iloc[0].round(5))
    Diff8.append(Av_NN['b9'].iloc[0].round(5))

    Diff9.append(Av_NS['a2'].iloc[0].round(5))
    Diff10.append(Av_NS['a3'].iloc[0].round(5))
    Diff11.append(Av_NS['a4'].iloc[0].round(5))
    Diff12.append(Av_NS['a5'].iloc[0].round(5))
    Diff13.append(Av_NS['a6'].iloc[0].round(5))
    Diff14.append(Av_NS['a7'].iloc[0].round(5))
    Diff15.append(Av_NS['a8'].iloc[0].round(5))
    Diff16.append(Av_NS['a9'].iloc[0].round(5))
    
    # print(abs(int(Av_NN['b2'].iloc[0])))
    mul=mul+0.01
    if Phi>2*np.pi:
         a=1
    if abs(int(Av_NN['b2'].iloc[0]))<1:
        
        if abs(int(Av_NN['b3'].iloc[0]))<1:
            
            if int(Av_NS['a2'].iloc[0])>2500:
                print('a2 pos',Phi/np.pi)
            elif int(Av_NN['a2'].iloc[0])<-2500:
                print('a2 neg',Phi/np.pi)
                    

plt.figure()
plt.plot(Angle,Diff1,label='b2')
plt.plot(Angle,Diff2,label='b3')
plt.plot(Angle,Diff3,label='b4')
plt.plot(Angle,Diff4,label='b5')
plt.plot(Angle,Diff5,label='b6')
plt.plot(Angle,Diff6,label='b7')
plt.plot(Angle,Diff7,label='b8')
plt.plot(Angle,Diff8,label='b9')
plt.hlines(0,0,2)
plt.title('Multipoles as a funtion of the axis rotation angle')
plt.xlabel('Angle (pi rad)')
plt.ylabel('bn')
plt.legend()    
    

plt.figure()
plt.plot(Angle,Diff9,label='a2')
plt.plot(Angle,Diff10,label='a3')
plt.plot(Angle,Diff11,label='a4')
plt.plot(Angle,Diff12,label='a5')
plt.plot(Angle,Diff13,label='a6')
plt.plot(Angle,Diff14,label='a7')
plt.plot(Angle,Diff15,label='a8')
plt.plot(Angle,Diff16,label='a9')

plt.hlines(0,0,2)
plt.title('Multipoles as a funtion of the axis rotation angle')
plt.xlabel('Angle (pi rad)')
plt.ylabel('an')
plt.legend() 








                



    






    
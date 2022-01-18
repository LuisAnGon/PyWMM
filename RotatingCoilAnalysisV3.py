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
import collections
from functionsV3 import RotatingCoilAnalysisTurn,Parameters, ContinuousRotatingCoilAnalysis, GetSensitivities,ReadFromC,Differences
from Postprocessing_Functions_RCAV3 import ReadRoxie, interpolate_Roxie_MM


step=65

#folder=r'C:\Users\Luis González\cernbox\Work\CIEMAT-\MCBXF\FFMM\rotating coil processing\Python\Data at 927\Scan-Fringe'
# folder=r'C:\Magnetic Measurements\Single Coils Centre\Up'

folder=r'C:\Magnetic Measurements\Training TE-MSC-MM\Scans\Single Dipoles Scans\MCBXFBP2-OD-test'
nombre=folder.split("\\")[-1]
print(folder)
# senspath=folder+'\\Kn_R45_130_N1_0001_A_AC.txt'   
senspath=r'Kn_R45_130_N1_0001_A_AC.txt'  
[knAbs,knCmp]=GetSensitivities(senspath)

meas=[] #Each folder Corresponds to one measurement in one position
dic={}        # We have to make sure that we are reading the positions in the correct order

separator="pos"

for name in os.listdir(folder):
    if os.path.isdir(os.path.join(folder, name)):
        if len(name.split(separator)[-1])==1:
            dic["0"+name.split(separator)[-1]]=name #in order to sort the files we need to avoid to have keys with only one digit
        else:
            dic[name.split(separator)[-1]]=name

keylist = list(dic.keys())
keylist.sort()


Av_pos=pd.DataFrame()
Av_neg=pd.DataFrame()

i_pos=1
paso=[]

for key in keylist:
#for measurement in meas:  # DIFFERENT MEASUREMENTS Meanining measurements at different positions
    #print(measurement)
    paso.append(step*(i_pos-1))
    
    for item in os.listdir(folder+'\\'+dic[key]): #Get the parameters from the exported file from C++
        if item.split('_')[-1]=='Parameters.txt':# Read the parameters at each position
            parameterpath=folder+'\\'+dic[key]+'\\'+item
            [p_turn, num_FDI,MagOrder,Rref,AnalysisOptions]=Parameters(parameterpath)
    
    
    for dire in [x[0] for x in os.walk(folder+'\\'+dic[key])][1:]:# We go to each position
        
        Current=dire.split('_')[-1] # Read the current at which the measurement was done
       
        for doc in os.listdir(dire): #We obtain every measurment done at each position
            #print('doc',doc)
            if 'raw_measurement_data' in doc:
                datafile=dire+'\\'+doc 
                
                Av=ContinuousRotatingCoilAnalysis (datafile, p_turn, knAbs, knCmp, MagOrder, Rref, AnalysisOptions,step) #Get the average of all the Rotation coil turns
                         
    
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
# Av_neg.loc["position"]=paso
#Av_neg.loc["B1"]=-Av_neg.loc["B1"]
# Av_pos.loc["B1"]=abs(Av_pos.loc["B1"]) #Make sure that the average of B1 is always positive)
# Av_pos.loc["position"]=paso

Av_neg.to_excel(folder+"\\neg.xlsx")
Av_pos.to_excel(folder+"\\pos.xlsx")


Av=pd.DataFrame()
for col in Av_neg.columns:
    Av[col]=(Av_neg[col]+Av_pos[col])/(2)
    
    if np.sign(Av_neg[col].loc["B1"]) == np.sign(Av_pos[col].loc["B1"]):
        Av[col].loc["B1"] = abs((Av_neg[col].loc["B1"] + Av_pos[col].loc["B1"])/2)
    else:
        Av[col].loc["B1"] = -(Av_pos[col].loc["B1"]-Av_neg[col].loc["B1"])/2
        
    
    # Av[col].loc["B1"]=(max(Av_neg[col].loc["B1"],Av_pos[col].loc["B1"])-min(Av_neg[col].loc["B1"],Av_pos[col].loc["B1"]))/2
    # print("max: ",max(Av_neg[col].loc["B1"],Av_pos[col].loc["B1"]))
    # print("min: ",min(Av_neg[col].loc["B1"],Av_pos[col].loc["B1"]))
    # print("Average: ",Av[col].loc["B1"])
    
Av.to_excel(folder+"\\Summary_MM_"+nombre+".xlsx")

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
                   
Av_NN['position']=Av_NN['position']-(Av_NN['position'].max()/2)
Av_NS['position']=Av_NS['position']-(Av_NS['position'].max()/2)


# ATTENTION MUST DO THIS FOR THE SKEW MAGNETS TOO!!!
# ATTENTION MUST DO THIS FOR THE SKEW MAGNETS TOO!!!
#Av_NN["B1"]=10000*Av_NN["B1"]/Av_NN["B1"].iloc[int(max(Av_NN.index)/2)]
# ATTENTION MUST DO THIS FOR THE SKEW MAGNETS TOO!!!
# ATTENTION MUST DO THIS FOR THE SKEW MAGNETS TOO!!!


#### READ OUTPUT FROM C++
# [summarypos,summaryneg]=ReadFromC(folder+'\\_Average_results.txt',step)


# # Find differences between Python and C++ Results

# dif_pos=Differences(Av_pos,summarypos)
# dif_neg=Differences(Av_neg,summaryneg)

# print(dif_pos)
# print(dif_neg)


[RoxieMP,Roxieall]=ReadRoxie(folder)

# print('ROXIE ALL:',RoxieMP)
if sum(RoxieMP.index.isin(['A16 Roxie']))>0:
    RoxieMP=RoxieMP.drop('A16 Roxie',axis=0)
    
RoxieMP.index=["B1","b2","b3","b4","b5","b6","b7","b8","b9","b10","b11","b12","b13","b14","b15",
               "A1","a2","a3","a4","a5","a6","a7","a8","a9","a10","a11","a12","a13","a14","a15"]
pos=Av.loc[["B1","b2","b3","b4","b5","b6","b7","b8","b9","b10","b11","b12","b13","b14","b15",
            "A1","a2","a3","a4","a5","a6","a7","a8","a9","a10","a11","a12","a13","a14","a15"]][Av.columns[int(len(Av.columns)/2)]]
   
    
comp_Roxie_Meas=pd.concat([RoxieMP,pos],axis=1)
comp_Roxie_Meas=comp_Roxie_Meas.round(6)
comp_Roxie_Meas.columns=['Roxie','MM']
print(comp_Roxie_Meas)

comp_Roxie_Meas.to_excel(folder+'\\Roxie_vs_MM_MiddlePoint '+nombre+'.xlsx')

roxie_int_Av_NN=interpolate_Roxie_MM(folder,Roxieall,Av_NN,'Av_NN')
roxie_int_Av_NS=interpolate_Roxie_MM(folder,Roxieall,Av_NS,'Av_NS')


# ATTENTION MUST DO THIS FOR THE SKEW MAGNETS TOO!!!
# ATTENTION MUST DO THIS FOR THE SKEW MAGNETS TOO!!!
for col in roxie_int_Av_NN.columns[2:16]:
    roxie_int_Av_NN[col]=10000*roxie_int_Av_NN[col]/roxie_int_Av_NN["B1 Roxie"]
roxie_int_Av_NN["B1 Roxie"]=10000*roxie_int_Av_NN["B1 Roxie"]/roxie_int_Av_NN["B1 Roxie"].iloc[int(max(roxie_int_Av_NN.index)/2)]
# ATTENTION MUST DO THIS FOR THE SKEW MAGNETS TOO!!!
# ATTENTION MUST DO THIS FOR THE SKEW MAGNETS TOO!!!


Roxie_MM_NN=pd.merge(Av_NN,roxie_int_Av_NN)
Roxie_MM_NN.to_excel(folder+"\\Roxie_vs_MM_Normal"+nombre+".xlsx")
Roxie_MM_NS=pd.merge(Av_NS,roxie_int_Av_NS)
Roxie_MM_NS.to_excel(folder+"\\Roxie_vs_MM_Skew"+nombre+".xlsx")

# plt.figure()
# plt.plot(Roxie_MM_NN["position"],Roxie_MM_NN["B1 Roxie"],marker="o",label="B1 Roxie Normalized to 0mm")
# plt.plot(Roxie_MM_NN["position"],Roxie_MM_NN["B1"],marker="o",label="B1 Normalized to 0mm")
# plt.legend()
# plt.xlabel("Position [mm]")
# plt.ylabel("Units")


# for num in range(2,16):
#     plt.figure()
#     plt.plot(Roxie_MM_NN["position"],Roxie_MM_NN["B"+str(num)+" Roxie"],marker="o",label="B"+str(num)+" Roxie")
#     plt.plot(Roxie_MM_NN["position"],Roxie_MM_NN["b"+str(num)],marker="o",label="b"+str(num))
#     plt.legend()
#     plt.xlabel("Position [mm]")
#     plt.ylabel("Units")
    
    # plt.figure()
    # plt.plot(Roxie_MM_NS["position"],Roxie_MM_NS["A"+str(num)+" Roxie"])
    # plt.plot(Roxie_MM_NS["position"],Roxie_MM_NS["a"+str(num)])




                



    






    
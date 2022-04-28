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
from Postprocessing_Functions_RCA import RotatingCoilAnalysisTurn, ContinuousRotatingCoilAnalysis, ReadRoxie, interpolate_Roxie_MM, GetSensitivities,SelectRoxie


# =============================================================================
# Takes manually the location of the folder with the data In each folder there are sub folders corresponding to each position along the magnet
# =============================================================================
#Folder used for measurements
#folder=r'C:\RCA\WMM_All_configurations\MCBXFB_02_Outer_Iron_Normal_600_20210510' 
#Folder used for developing
folder=r'C:\Users\admin\cernbox\Work\CIEMAT-\MCBXF\FFMM\Python Post-processing\PyWMM\WMM_All_configurations\MCBXFA_02_Outer_Collar_600_20220510'
nombre=folder.split("\\")[-1] #Folder name
print("Dipole measured: ", nombre.split("_")[0], nombre.split("_")[2], nombre.split("_")[3])
print("Rotating Coil Length: ", nombre.split("_")[-2])
# =============================================================================
# Takes manually the location of the folder with the data In each folder there are sub folders corresponding to each position along the magnet
# =============================================================================


# =============================================================================
# From the name of the folder defines: Length of Coil, Whether it is MCBXFA or MCBXFB
# =============================================================================
# AnalysisOptions= "cel deb dri rot"
AnalysisOptions= "cel deb dri nor rot"
MagOrder=1
Rref=0.05
p_turn=512
coil=nombre.split("_")[4]      #Detects the length of the coil from the name of the folder
tipo=nombre.split("_")[0][-1] # Detects wether Magnet is A or B

nombreShort=nombre.split("_200")[0].split("_600")[0]
iron=int(nombre.split("_")[3]=="Iron")
collar=int(nombre.split("_")[3]=="Collar")
inner=int(nombre.split("_")[2]=="Inner")
outer=int(nombre.split("_")[2]=="Outer")
# =============================================================================
# From the name of the folder defines: Length of Coil, Whether it is MCBXFA or MCBXFB
# =============================================================================


# =============================================================================
# Define the main orientation of the field based on wether it is inner or aouter / Normal or Skew
# =============================================================================


if collar==0 and iron==1:
    Rot90=0 
    print("The dipole is NOT rotated 90 deg")
if collar==1 and iron==0 and inner==1 and outer==0:
    Rot90=0
    print("The dipole is NOT rotated 90 deg")
if collar==1 and iron==0 and inner==0 and outer ==1:
    Rot90=1
    print("The dipole is rotated 90 deg")
    
if Rot90==1 and inner==1:
    MainField="SKEW"
elif Rot90==0 and inner==0:
    MainField="SKEW"
else:
    MainField="NORMAL"
# =============================================================================
# Define the main orientation of the field based on wether it is inner or aouter / Normal or Skew
# =============================================================================


print("Main Field: ",MainField)



# =============================================================================
# Based on above defines lists with the longitudinal positions (paso)
# =============================================================================
if coil=="600":
    if tipo == "A":
        paso=[0,600,1200,1700,2200,2800,3400]
        paso_middle=[1200,1700,2300]
    elif tipo == "B":
        paso=[0,600,1200,1800,2400]
        paso_middle=[1200]
elif coil == "200":
    if tipo == "A":
        paso=[200,  400,  600,  800, 1000, 1200, 1400, 1600, 1800, 2000,
        2200, 2400, 2600, 2800]
        paso_middle=[1000, 1200, 1400, 1600, 1800, 2000]
    elif tipo == "B":
        paso=[200,  400,  600,  800, 1000, 1200, 1400,1600,1800]
        paso_middle=[800, 1000, 1200]
# =============================================================================
# Based on above defines lists with the longitudinal positions (paso)
# =============================================================================


# =============================================================================
# Reads the sensitivities from The folder containig the .py code
# =============================================================================
senspath=r'C:\RCA\rca_calibration_data\Kn_R45_130_N1_0001_A_AC.txt'  
[knAbs,knCmp]=GetSensitivities(senspath)
# =============================================================================
# Reads the sensitivities from The folder containig the .py code
# =============================================================================



# =============================================================================
# From the data folder selects: The fluxes data .txt (For all positions along the magnet)   and The raw data (only one file)
# =============================================================================
meas=[] # Each folder Corresponds to one measurement in one value of the current for each position
dic={}  # We have to make sure that we are reading the positions in the correct order
Av_pos=pd.DataFrame()
Av_neg=pd.DataFrame()

i_pos=1

for name in os.listdir(folder):
   if name.split(".")[-1]=="txt":
       if name.split("_")[-1]=="Ascii.txt":
           meas.append(name)
       elif name.split(".")[0].split("_")[-1]=="raw":
           rawname=name
           rawfile=pd.read_csv(folder+"\\"+name, sep='\t')

raw=rawfile.groupby(["PosN","Iset(A)"]).first()
raw.to_excel(folder+"\\"+rawname+".xlsx")
# =============================================================================
# From the data folder selects: The fluxes data .txt (For all positions along the magnet)   and The raw data (only one file)
# =============================================================================





# =============================================================================
# Analysis of each position along the magnet:
# =============================================================================
for m in meas:# We go to each position
    pos=int(m.split("pos")[1].split("_")[0])
    Current=int(m.split("_")[9].split("A")[0]) # Read the current at which the measurement was done   
    step=int(m.split("_")[4].split("A")[0])
    
    file=folder+"\\"+m
    Av=ContinuousRotatingCoilAnalysis (file, p_turn, knAbs, knCmp, MagOrder, Rref, AnalysisOptions,step,MainField) #Get the average of all the Rotation coil turns
    
# =============================================================================
#  Sets the positions “paso” as the column “position” and Concatenates df to Av_pos or Av_neg depending on the signal of the current
# =============================================================================
    Av=Av.rename(columns={"Average":i_pos-1})
    
    if Current>0:
        Av_pos=pd.concat([Av_pos,Av],axis=1)
        
        
    if Current<0: 
        Av_neg=pd.concat([Av_neg,Av],axis=1)
        
# =============================================================================
#  Sets the positions “paso” as the column “position” and Concatenates df to Av_pos or Av_neg depending on the signal of the current
# =============================================================================
        





# =============================================================================
# Analysis of each position along the magnet:
# =============================================================================

# =============================================================================
# Creates Av which is the average of the results for Positive and Negative current and saves it as Summary.
# =============================================================================
Av_neg.columns=paso    
# Specify the position of each measurement according to the number of position
Av_neg.loc["position"]=paso

Av_pos.columns=paso   
Av_pos.loc["position"]=paso

Av_neg.astype(float)
Av_pos.astype(float)


Av=pd.DataFrame()

if MainField=="NORMAL":
    mainRoxie="B1"
elif MainField=="SKEW":
    mainRoxie="A1"

for col in Av_neg.columns:
    
    if np.sign(Av_neg[col].loc[mainRoxie]) == 1: 
        Av_neg[col].loc[mainRoxie]=-Av_neg[col].loc[mainRoxie]
    elif np.sign(Av_pos[col].loc[mainRoxie]) == -1:
        Av_neg[col].loc[mainRoxie]=-Av_neg[col].loc[mainRoxie]
        
    Av[col]=(Av_pos[col]+Av_neg[col])/(2)
# =============================================================================
# Creates Av which is the average of the results for Positive and Negative current and saves it as Summary.    
# =============================================================================





# =============================================================================
# Takes the A´s, a´s, B´s and b´s and separates Av in Av Normal and Av Skew
# =============================================================================
Av_NN=pd.DataFrame()
Av_NS=pd.DataFrame()



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
# =============================================================================
# Takes the A´s, a´s, B´s and b´s and separates Av in Av Normal and Av Skew
# =============================================================================





# =============================================================================
# Centers the position of the measurements to be from negative to positive
# =============================================================================
if MainField=="SKEW":
    Av_NS[mainRoxie]=10000*Av_NS[mainRoxie]/Av_NS[mainRoxie][(Av_NS.index.max()-Av_NS.index.min())/2]
elif MainField=="NORMAL":
    Av_NN[mainRoxie]=10000*Av_NN[mainRoxie]/Av_NN[mainRoxie][(Av_NN.index.max()-Av_NN.index.min())/2]
# =============================================================================
# Centers the position of the measurements to be from negative to positive
# =============================================================================


# Av.rename(index={mainRoxie:"Main Field [mT]"},inplace=True)
Av.loc["Current [A]"]=[5]*len(paso)
Av.loc["Main Field [mT]"]=Av.loc[mainRoxie]*1000
Av.loc["Transfer Function [mT/A"]=Av.loc["Main Field [mT]"]/Av.loc["Current [A]"]



integrField=sum(Av.loc["Main Field [mT]"])*(int(coil)/1000)
Av.loc["%"]=Av.loc["Main Field [mT]"]*int(coil)/1000/integrField


Av.to_excel(folder+"\\Summary_MM_"+nombre+".xlsx")




#*****************************************************************************************************
#****************************************COMPARISON WITH ROXIE****************************************
#*****************************************************************************************************


# =============================================================================
# Finds the Roxie file in C: according to the characteristics of the magnet. Info taken from folder name
# =============================================================================
RoxieFolder=SelectRoxie(iron,inner,Rot90)
print("Roxie Folder: ", RoxieFolder)
# =============================================================================
# Finds the Roxie file in C: according to the characyeristics of the magnet. Info taken from folder name
# =============================================================================



# =============================================================================
# Reads the simulations performed with Roxie
# =============================================================================
[RoxieMP,Roxieall]=ReadRoxie(RoxieFolder,norm=False,skew=int(MainField=="SKEW"))

#If We have a long magnet MCBXFA, the Roxie file must be extended +1m
if tipo=="A":
    newposition=[]
    positionneg=list(Roxieall.position[Roxieall.position<0]-500)
    positionpos=list(Roxieall.position[Roxieall.position>0]+500)
    newposition=positionneg+positionpos
    Roxieall.position=newposition
    




# =============================================================================
# Reads the simulations performed with Roxie
# =============================================================================



# =============================================================================
# Removes the path along the field in case it was calculated in Roxie
# =============================================================================
if sum(RoxieMP.index.isin(['A16 Roxie']))>0:
    RoxieMP=RoxieMP.drop('A16 Roxie',axis=0)
# =============================================================================
# Removes the path along the field in case it was calculated in Roxie
# =============================================================================
    


# =============================================================================
# With MP Creates Comp_Roxie_Meas: A df to compare the values of Roxie in the middle and the measurement in the middle. Then saves it
# =============================================================================
RoxieMP.index=["B1","b2","b3","b4","b5","b6","b7","b8","b9","b10","b11","b12","b13","b14","b15",
               "A1","a2","a3","a4","a5","a6","a7","a8","a9","a10","a11","a12","a13","a14","a15"]
pos=Av.loc[["B1","b2","b3","b4","b5","b6","b7","b8","b9","b10","b11","b12","b13","b14","b15",
            "A1","a2","a3","a4","a5","a6","a7","a8","a9","a10","a11","a12","a13","a14","a15"]][Av.columns[int(len(Av.columns)/2)]]
   
    
comp_Roxie_Meas=pd.concat([RoxieMP,pos],axis=1)
comp_Roxie_Meas=comp_Roxie_Meas.round(6)
comp_Roxie_Meas.columns=['Roxie','MM']
print("Comp_Roxie_Meas",comp_Roxie_Meas)

comp_Roxie_Meas.to_excel(folder+'\\Roxie_vs_MM_MiddlePoint '+nombre+'.xlsx')
# =============================================================================
# With MP Creates Comp_Roxie_Meas: A df to compare the values of Roxie in the middle and the measurement in the middle. Then saves it
# =============================================================================



# =============================================================================
# Creates a df in which the Roxie profiles are interpolated to select the values corresponding to the measured points and is merged to the measured points. 
# Then Saves a different file for normal and skew 
# =============================================================================

if mainRoxie=="B1":
    roxie_int_Av_NN=interpolate_Roxie_MM(folder,Roxieall,Av_NN,'Av_NN',tipo)
    roxie_int_Av_NN=roxie_int_Av_NN[["position","B1 Roxie","B2 Roxie","B3 Roxie","B4 Roxie","B5 Roxie","B6 Roxie","B7 Roxie","B8 Roxie",
                                      "B9 Roxie","B10 Roxie","B11 Roxie","B12 Roxie","B13 Roxie","B14 Roxie","B15 Roxie"]]
    Roxie_MM_NN=pd.merge(Av_NN,roxie_int_Av_NN)
    
    
    for b in ['b2', 'b3', 'b4', 'b5', 'b6', 'b7', 'b8', 'b9', 'b10','b11', 'b12', 'b13', 'b14', 'b15']:
        Roxie_MM_NN[b]=Roxie_MM_NN[b]*Roxie_MM_NN["B1"]/10000
    Roxie_MM_NN.to_excel(folder+"\\Roxie_vs_MM_Normal"+nombre+".xlsx")


elif mainRoxie=="A1":
    roxie_int_Av_NS=interpolate_Roxie_MM(folder,Roxieall,Av_NS,'Av_NS',tipo)
    roxie_int_Av_NS=roxie_int_Av_NS[["position","A1 Roxie","A2 Roxie","A3 Roxie","A4 Roxie","A5 Roxie","A6 Roxie","A7 Roxie","A8 Roxie",
                                      "A9 Roxie","A10 Roxie","A11 Roxie","A12 Roxie","A13 Roxie","A14 Roxie","A15 Roxie"]]
    Roxie_MM_NS=pd.merge(Av_NS,roxie_int_Av_NS)
    
    for b in ['b2', 'b3', 'b4', 'b5', 'b6', 'b7', 'b8', 'b9', 'b10','b11', 'b12', 'b13', 'b14', 'b15']:
        Roxie_MM_NN[b]=Roxie_MM_NN[b]*Roxie_MM_NN["B1"]/10000
    
    Roxie_MM_NS.to_excel(folder+"\\Roxie_vs_MM_Skew"+nombre+".xlsx")

# =============================================================================
# Creates a df in which the Roxie profiles are interpolated to select the values corresponding to the measured points and is merged to the measured points. 
# Then Saves a different file for normal and skew 
# =============================================================================



#*****************************************************************************************************
#****************************************PLOTTING****************************************
#*****************************************************************************************************


fig, axs = plt.subplots(2,figsize=(5,7))
x=[2,3,4,5,6,7,8,9,10,11,12,13,14,15]
axs[0].bar([a-0.15 for a in x],comp_Roxie_Meas.loc[['b2', 'b3', 'b4', 'b5', 'b6', 'b7', 'b8', 'b9', 'b10', 'b11','b12', 'b13', 'b14', 'b15']]["Roxie"],label="Roxie",width=0.3)
axs[0].bar([a+0.15 for a in x],comp_Roxie_Meas.loc[['b2', 'b3', 'b4', 'b5', 'b6', 'b7', 'b8', 'b9', 'b10', 'b11','b12', 'b13', 'b14', 'b15']]["MM"],label="MM",width=0.3)
axs[0].set_title("Comparison Roxie-MM "+nombreShort)
axs[0].grid(linestyle='--', linewidth=0.5)
axs[0].set_xticks([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16])
axs[0].set_yticks([-20,-15,-10,-5,0,5,10,15,20])
#axs[0].set_xlabel("n")
axs[0].set_ylabel("bn")
axs[0].legend()


axs[1].bar([a-0.15 for a in x],comp_Roxie_Meas.loc[['a2', 'a3', 'a4', 'a5', 'a6', 'a7', 'a8', 'a9', 'a10', 'a11','a12', 'a13', 'a14', 'a15']]["Roxie"],label="Roxie",width=0.3)
axs[1].bar([a+0.15 for a in x],comp_Roxie_Meas.loc[['a2', 'a3', 'a4', 'a5', 'a6', 'a7', 'a8', 'a9', 'a10', 'a11','a12', 'a13', 'a14', 'a15']]["MM"],label="MM",width=0.3)
#axs[1].set_title("Comparison Roxie-MM Skew Multipoles")
axs[1].grid(linestyle='--', linewidth=0.5)
axs[1].set_xticks([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16])
axs[1].set_yticks([-20,-15,-10,-5,0,5,10,15,20])
axs[1].set_xlabel("n")
axs[1].set_ylabel("an")
axs[1].legend()

fig.savefig((folder+"\\Informe_Multipoles_Center.pdf"))

if mainRoxie=="B1":
    u1="B1"
    u3="B3"
    u5="B5"
    l1="b1"
    l3="b3"
    l5="b5"
if mainRoxie=="A1":
    u1="A1"
    u3="A3"
    u5="A5"
    l1="a1"
    l3="a3"
    l5="a5"
    

fig, axs = plt.subplots(3,figsize=(5,7))
axs[0].plot(roxie_int_Av_NN.position,Roxie_MM_NN[u1+" Roxie"],label="Roxie "+u1,marker=".")
axs[0].plot(roxie_int_Av_NN.position,Roxie_MM_NN[u1],label="MM "+u1,marker=".")
axs[0].set_xticks(list(roxie_int_Av_NN.position))
axs[0].set_xticks([it+300 for it in list(roxie_int_Av_NN.position)]+[it-300 for it in list(roxie_int_Av_NN.position)],minor=True)
axs[0].xaxis.grid(False, which='major')
axs[0].xaxis.grid(True, linestyle='--', linewidth=0.5, which='minor')
axs[0].set_ylabel(u1)
axs[0].legend()

axs[1].plot(roxie_int_Av_NN.position,Roxie_MM_NN[u3+" Roxie"],label="Roxie "+l3,marker=".")
axs[1].plot(roxie_int_Av_NN.position,Roxie_MM_NN[l3],label="MM "+l3,marker=".")
axs[1].set_xticks(list(roxie_int_Av_NN.position))
axs[1].set_xticks([it+300 for it in list(roxie_int_Av_NN.position)]+[it-300 for it in list(roxie_int_Av_NN.position)],minor=True)
axs[1].xaxis.grid(False, which='major')
axs[1].xaxis.grid(True, linestyle='--', linewidth=0.5, which='minor')
axs[1].set_ylabel(l3)
axs[1].legend()

axs[2].plot(roxie_int_Av_NN.position,Roxie_MM_NN[u5+" Roxie"],label="Roxie "+l5,marker=".")
axs[2].plot(roxie_int_Av_NN.position,Roxie_MM_NN[l5],label="MM "+l5,marker=".")
axs[2].set_xticks(list(roxie_int_Av_NN.position))
axs[2].set_xticks([it+300 for it in list(roxie_int_Av_NN.position)]+[it-300 for it in list(roxie_int_Av_NN.position)],minor=True)
axs[2].xaxis.grid(False, which='major')
axs[2].xaxis.grid(True, linestyle='--', linewidth=0.5, which='minor')
axs[2].set_ylabel(l5)
axs[2].set_xlabel("Position [mm]")
axs[2].legend()



fig.savefig((folder+"\\Informe_Multipoles_Profile.pdf"))

# plt.figure()
# plt.plot(Roxie_MM_NN["position"],Roxie_MM_NN["B1 Roxie"],marker="o",label="B1 Roxie Normalized to 0mm")
# plt.plot(Roxie_MM_NN["position"],Roxie_MM_NN["B1"],marker="o",label="B1 Normalized to 0mm")
# plt.legend()
# plt.xlabel("Position [mm]")
# plt.ylabel("Units")


# for num in range(2,16):
    # plt.figure()
    # plt.plot(Roxie_MM_NN["position"],Roxie_MM_NN["B"+str(num)+" Roxie"],marker="o",label="B"+str(num)+" Roxie")
    # plt.plot(Roxie_MM_NN["position"],Roxie_MM_NN["b"+str(num)],marker="o",label="b"+str(num))
    # plt.legend()
    # plt.xlabel("Position [mm]")
    # plt.ylabel("Units")
    
    # plt.figure()
    # plt.plot(Roxie_MM_NS["position"],Roxie_MM_NS["A"+str(num)+" Roxie"])
    # plt.plot(Roxie_MM_NS["position"],Roxie_MM_NS["a"+str(num)])




                



    






    
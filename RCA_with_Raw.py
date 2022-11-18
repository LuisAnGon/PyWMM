# -*- coding: utf-8 -*-
"""
Created on Wed May 18 09:47:43 2022

@author: luis gonzalez
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
from Postprocessing_Functions_RCA_Raw import RotatingCoilAnalysisTurn, ContinuousRotatingCoilAnalysisRaw, ReadRoxie, interpolate_Roxie_MM, GetSensitivities,SelectRoxie,ReadShim


# =============================================================================
# Takes manually the location of the folder with the data. In the folder there is a raw.txt file and a shimming.output file

# folder=r'C:\PyWMM-Medidas Elytt\MCBXFB02\Outer Collar'
folder=r'C:\PyWMM-Medidas Elytt\MCBXFB02\Inner Iron'
nombre=folder.split("\\")[-1] #Folder nameprint
shimpath=r"C:\PyWMM\Shimming"


# Takes manually the location of the folder with the data. In the folder there is a raw.txt file and a shimming.output file
# =============================================================================

for name in os.listdir(folder):
   if name.split("_")[-1]=="raw.txt": #Selects the raw file
       rawname=folder+"\\"+name
   elif name.split(".")[-1]=="output":#Selects the shimming file
       shimname=folder+"\\"+name


# =============================================================================
#Reads the raw file and converts it into a dataframe. A bit tricky since in the .txt there are columns without names 
#(all those corresponding to fluxes ans time)
newrawname=folder+"\\newraw.txt"

f=open(rawname)
print("rawname: ",rawname)
head=f.readlines()[0]
f.close()
head=head.split("\t")



# delete the header line
with open(rawname, "r") as f:
    lines = f.readlines()
with open(newrawname, "w") as ff:
    for line in lines[1:]:
        ff.write(line)

rawfile=pd.read_csv(newrawname, sep='\t', header=None)

os.remove(newrawname)

headerlst=list(rawfile.columns)

newheader=head+headerlst[len(head):]

# add the new header line
rawfile.columns=newheader

rawfile.to_excel(folder+"\\Excel_Raw.xlsx")

#Reads the raw file and converts it into a dataframe. A bit tricky since in the .txt there are columns without names 
#(all those corresponding to flufes ans time)
# =============================================================================



# =============================================================================
# From the Raw file defines parameters

AnalysisOptions= "cel deb dri nor rot fed"
MagOrder=1
Rref=0.05
p_turn=512
CalibAng=23.2
coil=str(rawfile['Lcoil(m)'][0])     #Detects the length of the coil from the name of the folder
print("coil",coil)
tipo=rawfile['Magnet type'][0][-1] # Detects wether Magnet is A or B
print("tipo",tipo)

# nombreShort=nombre.split("_200")[0].split("_600")[0]
iron=int(rawfile['Fabrication stage'][0]=="Iron")
print('Fabrication stage',rawfile['Fabrication stage'][0])
collar=int(rawfile['Fabrication stage'][0]=="Collar")
inner=int(rawfile['Dipole'][0]=="Inner")
outer=int(rawfile['Dipole'][0]=="Outer")
print("inner-outer",rawfile['Dipole'][0])

# From the Raw file defines parameters
# =============================================================================


# =============================================================================
# Define the main orientation of the field based on wether it is inner or outer / Normal or Skew

if collar==0 and iron==1:
    Rot90=0 
    print("The dipole is NOT rotated 90 deg")
if collar==1 and iron==0 and inner==1 and outer==0:
    Rot90=0
    print("The dipole is NOT rotated 90 deg")
if collar==1 and iron==0 and inner==0 and outer ==1:
    Rot90=1
    print("The dipole is rotated 90 deg")

Rot90=int(input("Rot90"))
    
if Rot90==1 and inner==1:
    MainField="SKEW"
elif Rot90==0 and inner==0:
    MainField="SKEW"
else:
    MainField="NORMAL"

# MainField=str(input("Mainfield"))

# Define the main orientation of the field based on wether it is inner or outer / Normal or Skew
# =============================================================================

print("coil",coil)
print("Main Field: ",MainField)

# =============================================================================
# Based on above defines lists with the longitudinal positions (paso)y

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

# Based on above defines lists with the longitudinal positions (paso)
# =============================================================================



# =============================================================================
# Reads the sensitivities from The folder containig the .py code

senspath=r'C:\PyWMM\rca_calibration_data\Kn_DQS_5_24_16_250_115x650_0001_A_AC.txt'  
[knAbs,knCmp]=GetSensitivities(senspath)

# Reads the sensitivities from The folder containig the .py code
# =============================================================================


# =============================================================================
# Goes through the raw file and makes a continuiuous
# =============================================================================

poscurr=[]#list of positive current values. Updated for each position and rotation direction
negcurr=[]#list of positive current values. Updated for each position and rotation direction
Av_pos=pd.DataFrame()
Av_neg=pd.DataFrame()
Av_pos_rot=pd.DataFrame()
Av_neg_rot=pd.DataFrame()
i_pos=1




for pos in list(rawfile['PosN'].unique()):
    df=rawfile[rawfile["PosN"]==pos]
    # print("pos:" ,pos)
    for Current in list(rawfile['Iset(A)'].unique()):
        dfc=df[df["Iset(A)"]==Current]
        # print("Current: " ,Current)
        for rot in list(rawfile['Speed(rpm)'].unique()):
            # print("Rotation: " ,rot)
            dfr=dfc[dfc["Speed(rpm)"]==rot]
            # print(dfr["Imeas(A)"])
        
            Av=ContinuousRotatingCoilAnalysisRaw (dfr, knAbs, knCmp, MagOrder, Rref, AnalysisOptions,coil,MainField,rot,CalibAng) #Get the average of all the 10 Rotation coil turns for this pos,curr y rot_Direction
            # print (Av)
        # =============================================================================
        #  Concatenates df to Av_pos or Av_neg depending on the signal of the current
        
            Av=Av.rename(columns={"Average":i_pos-1})
            
            
    
            if Current>0:
                Av_pos=pd.concat([Av_pos,Av],axis=1) #For positive currents concatenates all turns at different positions and rotation directions
                # print(Av_pos)
                poscurr=poscurr+[float(dfr["Imeas(A)"].mean())] #Gets, FROM dfr, the average of the currents for all the 10 Rotation coil turns for this pos,curr y rot_Direction and concatenates in a list
                # print(poscurr)
                    
                    
            if Current<0: 
                Av_neg=pd.concat([Av_neg,Av],axis=1)#For negative currents concatenates all turns at different positions and rotation directions
                negcurr=negcurr+[-1*(float(dfr["Imeas(A)"].mean()))]  #Gets, FROM dfr, the average of the currents for all the 10 Rotation coil turns for this pos,curr y rot_Direction and concatenates in a list
            
            
    i_pos+=1
        
        #  Concatenates df to Av_pos or Av_neg depending on the signal of the current
        # =============================================================================



#Average          
Av_pos.loc["Current [A]"]=poscurr #adds the row current to the Av (For positive currents)
i=0

#Av_pos and Av_neg have this shape:
    #
    # |    Pos 0    |    Pos 1    |    Pos 2    |    Pos 3    |    Pos 4    |
    # |rot-60|rot 60|rot-60|rot 60|rot-60|rot 60|rot-60|rot 60|rot-60|rot 60|
    #


for item in paso: #for loop in steps of two in order to step from position to position. 
                  #Adds to the emptyAv_pos_rot df columns for each position and selects the multipoles for positive an negative rotations and makes the average
    Av_pos_rot[item]=(Av_pos.iloc[:,i]+Av_pos.iloc[:,i+1])/2
    i=i+2

Av_neg.loc["Current [A]"]=negcurr #adds the row current to the Av (For negative currents)
i=0

for item in paso: #for loop in steps of two in order to step from position to position.
                  #Adds to the emptyAv_pos_rot df columns for each position and selects the multipoles for positive an negative rotations and makes the average
    Av_neg_rot[item]=(Av_neg.iloc[:,i]+Av_neg.iloc[:,i+1])/2
    i=i+2


# =============================================================================
# Analysis of each position along the magnet:
# =============================================================================

# =============================================================================
# Creates Av which is the average of the results for Positive and Negative current and saves it as Summary.
# =============================================================================
# Av_neg.columns=paso    
# Specify the position of each measurement according to the number of position
Av_neg_rot.loc["position"]=paso

# Av_pos.columns=paso   
Av_pos_rot.loc["position"]=paso

Av_neg_rot.astype(float)
Av_pos_rot.astype(float)


Av=pd.DataFrame()

if MainField=="NORMAL":
    mainRoxie="B1"
elif MainField=="SKEW":
    mainRoxie="A1"

for col in Av_neg_rot.columns:
    
    if np.sign(Av_neg_rot[col].loc[mainRoxie]) == 1: 
        Av_neg_rot[col].loc[mainRoxie]=-Av_neg_rot[col].loc[mainRoxie]
    elif np.sign(Av_pos_rot[col].loc[mainRoxie]) == -1:
        Av_neg_rot[col].loc[mainRoxie]=-Av_neg_rot[col].loc[mainRoxie]
        
    Av[col]=(Av_pos_rot[col]+Av_neg_rot[col])/(2)
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
# Av.loc["Current [A]"]=[5]*len(paso)
Av.loc["Main Field [mT]"]=Av.loc[mainRoxie]*1000
Av.loc["Transfer Function [mT/A]"]=Av.loc["Main Field [mT]"]/Av.loc["Current [A]"]





integrField=sum(Av.loc["Main Field [mT]"])*(int(coil)/1000)
Av.loc["%"]=Av.loc["Main Field [mT]"]*int(coil)/1000/integrField

#Last correction of the angle due to the fact that the FFMM data has a cyclic shift of +1 so we must correct the effect of such +1 value of the encoder, which is (2*np.pi)/512
Angle_Corr=1000*(2*np.pi)/512
Av.loc["Angle (Applied Calib) [mrad]"]=Av.loc["Angle (Applied Calib) [mrad]"]-Angle_Corr

Av["Integral"]=[""]*len(Av)

for row in ["Current [A]","Main Field [mT]","b2","b3","b4","b5","b6","b7","b8","b9","b10","b11","b12","b13","b14","b15","a2","a3","a4","a5","a6","a7","a8","a9","a10","a11","a12","a13","a14","a15"]:
    if row=="Main Field [mT]":    
        Av["Integral"].loc[row]=sum(Av[paso[1:][:-1]].loc[row])*int(coil)/1000
    else:
        sp1=Av[paso[1:][:-1]].loc[row]
        sp2=Av[paso[1:][:-1]].loc["%"]
        Av["Integral"].loc[row]=np.dot(sp1,sp2)

Av["Integral"].loc["Transfer Function [mT/A]"]=Av["Integral"].loc["Main Field [mT]"]/Av["Integral"].loc["Current [A]"]
  

Av.to_excel(folder+"\\Summary_MM_"+nombre+".xlsx")




#*****************************************************************************************************
#****************************************COMPARISON WITH ROXIE****************************************
#*****************************************************************************************************


# =============================================================================
# Finds the Roxie file in C: according to the characteristics of the magnet. Info taken from folder name
# =============================================================================
RoxieFolder=SelectRoxie(iron,inner,Rot90)
print("Roxie Folder: ", RoxieFolder)

# RoxieFolder=r"C:\PyWMM\ROXIE-MCBXFB\No Iron\No Iron Inner Dipole\MCBXFB No Iron Inner Dipole Rot90"
# =============================================================================
# Finds the Roxie file in C: according to the characyeristics of the magnet. Info taken from folder name
# =============================================================================



# =============================================================================
# Reads the simulations performed with Roxie
# =============================================================================
[RoxieMP,Roxieall,RoxieInteg,RoxieMainInteg]=ReadRoxie(RoxieFolder,norm=False,skew=int(MainField=="SKEW"))

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
            "A1","a2","a3","a4","a5","a6","a7","a8","a9","a10","a11","a12","a13","a14","a15"]][Av.columns[int((len(Av.columns)/2)-1)]]
   
    
comp_Roxie_Meas=pd.concat([RoxieMP,pos],axis=1)
comp_Roxie_Meas=comp_Roxie_Meas.round(6)
comp_Roxie_Meas.columns=['Roxie_Middle','MM_Middle']

# print("Comp_Roxie_Meas:")
# print(comp_Roxie_Meas)

# comp_Roxie_Meas.to_excel(folder+'\\Roxie_vs_MM_MiddlePoint '+nombre+'.xlsx')
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
    
    for a in ['a2', 'a3', 'a4', 'a5', 'a6', 'a7', 'a8', 'a9', 'a10','a11', 'a12', 'a13', 'a14', 'a15']:
        Roxie_MM_NS[a]=Roxie_MM_NS[a]*Roxie_MM_NS["A1"]/10000
    
    Roxie_MM_NS.to_excel(folder+"\\Roxie_vs_MM_Skew"+nombre+".xlsx")

# =============================================================================
# Creates a df in which the Roxie profiles are interpolated to select the values corresponding to the measured points and is merged to the measured points. 
# Then Saves a different file for normal and skew 
# =============================================================================




#*****************************************************************************************************
#****************************************COMPARISON WITH SHIMMING (It is also ROXIE) ****************************************
#*****************************************************************************************************



# shimFolder=r"C:\PyWMM-Medidas Elytt\MCBXFB02\Inner Collar"

shimFolder=folder

[shimMP,shimall,shiminteg,shimmaininteg]=ReadRoxie(shimFolder,norm=False,skew=int(MainField=="SKEW"),bidi=True)

indexlist=["B1","b2","b3","b4","b5","b6","b7","b8","b9","b10","b11","b12","b13","b14","b15",
                "A1","a2","a3","a4","a5","a6","a7","a8","a9","a10","a11","a12","a13","a14","a15"]

# comp_Roxie_Meas_shim=pd.concat([comp_Roxie_Meas,shimMP],axis=1)
# comp_Roxie_Meas_shim.columns=["Roxie","MM","Shim"]

comp_Roxie_Meas_shim=comp_Roxie_Meas

comp_Roxie_Meas_shim["Shim_Middle"]=shiminteg

comp_Roxie_Meas_shim["Integrated_Roxie"]=RoxieInteg
comp_Roxie_Meas_shim["Integrated_MM"]=list(Av["Integral"].loc[indexlist])


if MainField=="NORMAL":
    comp_Roxie_Meas_shim["Integrated_Roxie"].loc["B1"]=RoxieMainInteg
    comp_Roxie_Meas_shim["Integrated_MM"].loc["B1"]=Av["Integral"].loc["Main Field [mT]"]
elif MainField=="SKEW":
    comp_Roxie_Meas_shim["Integrated_Roxie"].loc["A1"]=RoxieMainInteg
    comp_Roxie_Meas_shim["Integrated_MM"].loc["A1"]=Av["Integral"].loc["Main Field [mT]"]
    


# print("Comp_Roxie_Meas:")
# print(comp_Roxie_Meas_shim)
comp_Roxie_Meas_shim.to_excel(folder+'\\Roxie_vs_MM_vs_Shim_MiddlePoint '+nombre+'.xlsx')




#*****************************************************************************************************
#****************************************PLOTTING****************************************
#*****************************************************************************************************


fig, axs = plt.subplots(2,figsize=(5,7))
x=[2,3,4,5,6,7,8,9,10,11,12,13,14,15]
axs[0].bar([a-0.25 for a in x],comp_Roxie_Meas_shim.loc[['b2', 'b3', 'b4', 'b5', 'b6', 'b7', 'b8', 'b9', 'b10', 'b11','b12', 'b13', 'b14', 'b15']]["Shim_Middle"],label="Shimming",width=0.3,color="b")
axs[0].bar([a+0.25 for a in x],comp_Roxie_Meas_shim.loc[['b2', 'b3', 'b4', 'b5', 'b6', 'b7', 'b8', 'b9', 'b10', 'b11','b12', 'b13', 'b14', 'b15']]["Roxie_Middle"],label="Roxie",width=0.3,color="g")
axs[0].bar([a for a in x],comp_Roxie_Meas_shim.loc[['b2', 'b3', 'b4', 'b5', 'b6', 'b7', 'b8', 'b9', 'b10', 'b11','b12', 'b13', 'b14', 'b15']]["MM_Middle"],label="MM",width=0.5,color="r")

axs[0].set_title("Comparison Roxie-MM-Shimming ")
axs[0].grid(linestyle='--', linewidth=0.5)
axs[0].set_xticks([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16])
axs[0].set_yticks([-20,-15,-10,-5,0,5,10,15,20])
#axs[0].set_xlabel("n")
axs[0].set_ylabel("bn")
axs[0].legend()


axs[1].bar([a-0.25 for a in x],comp_Roxie_Meas_shim.loc[['a2', 'a3', 'a4', 'a5', 'a6', 'a7', 'a8', 'a9', 'a10', 'a11','a12', 'a13', 'a14', 'a15']]["Shim_Middle"],label="Shimming",width=0.3,color="b")
axs[1].bar([a+0.25 for a in x],comp_Roxie_Meas_shim.loc[['a2', 'a3', 'a4', 'a5', 'a6', 'a7', 'a8', 'a9', 'a10', 'a11','a12', 'a13', 'a14', 'a15']]["Roxie_Middle"],label="Roxie",width=0.3,color="g")
axs[1].bar([a for a in x],comp_Roxie_Meas_shim.loc[['a2', 'a3', 'a4', 'a5', 'a6', 'a7', 'a8', 'a9', 'a10', 'a11','a12', 'a13', 'a14', 'a15']]["MM_Middle"],label="MM",width=0.5,color="r")

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
    
if mainRoxie=="B1":
    roxie_int_Av_N=roxie_int_Av_NN
    Roxie_MM_N=Roxie_MM_NN
if mainRoxie=="A1":
    roxie_int_Av_N=roxie_int_Av_NS
    Roxie_MM_N=Roxie_MM_NS

fig, axs = plt.subplots(3,figsize=(5,7))
axs[0].plot(roxie_int_Av_N.position,Roxie_MM_N[u1+" Roxie"],label="Roxie "+u1,marker=".")
axs[0].plot(roxie_int_Av_N.position,Roxie_MM_N[u1],label="MM "+u1,marker=".")
axs[0].set_xticks(list(roxie_int_Av_N.position))
axs[0].set_xticks([it+300 for it in list(roxie_int_Av_N.position)]+[it-300 for it in list(roxie_int_Av_N.position)],minor=True)
axs[0].xaxis.grid(False, which='major')
axs[0].xaxis.grid(True, linestyle='--', linewidth=0.5, which='minor')
axs[0].set_ylabel(u1)
axs[0].legend()

axs[1].plot(roxie_int_Av_N.position,Roxie_MM_N[u3+" Roxie"],label="Roxie "+l3,marker=".")
axs[1].plot(roxie_int_Av_N.position,Roxie_MM_N[l3],label="MM "+l3,marker=".")
axs[1].set_xticks(list(roxie_int_Av_N.position))
axs[1].set_xticks([it+300 for it in list(roxie_int_Av_N.position)]+[it-300 for it in list(roxie_int_Av_N.position)],minor=True)
axs[1].xaxis.grid(False, which='major')
axs[1].xaxis.grid(True, linestyle='--', linewidth=0.5, which='minor')
axs[1].set_ylabel(l3)
axs[1].legend()

axs[2].plot(roxie_int_Av_N.position,Roxie_MM_N[u5+" Roxie"],label="Roxie "+l5,marker=".")
axs[2].plot(roxie_int_Av_N.position,Roxie_MM_N[l5],label="MM "+l5,marker=".")
axs[2].set_xticks(list(roxie_int_Av_N.position))
axs[2].set_xticks([it+300 for it in list(roxie_int_Av_N.position)]+[it-300 for it in list(roxie_int_Av_N.position)],minor=True)
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


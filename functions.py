# -*- coding: utf-8 -*-
"""
Created on Fri Apr 30 14:37:57 2021

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
def RotatingCoilAnalysisTurn(Fabs, Fcmp, knAbs, knCmp, MagOrder, Rref, AnalysisOptions):
    #fft
    f_abs=2*fft(np.array(Fabs))/len(Fabs)
    f_cmp=2*fft(np.array(Fcmp))/len(Fcmp)
    f_abs=f_abs[1:16]
    f_cmp=f_cmp[1:16]

    
    
    #Apply coil sensitivity
    Rref_vec=np.zeros(len(f_abs))
    
    for n in np.arange(0,len(Rref_vec)):
        Rref_vec[n]=Rref**n
    
    c_sens_abs = Rref_vec*(1/knAbs)*f_abs
    c_sens_cmp = Rref_vec*(1/knCmp)*f_cmp
    
    #Apply Rotation------------------------------------------------------------
    #Phase between +pi/2 and -pi/2
    
    SignalPhase=np.angle(c_sens_abs[(MagOrder-1)])
    if (SignalPhase>np.pi/2):
        SignalPhase = (SignalPhase-np.pi)
    elif (SignalPhase<-np.pi/2):
        SignalPhase = (SignalPhase+np.pi)
    
    #direction of the main component
    
    PhiOut=SignalPhase/MagOrder
    
    
    B_Main_Rotated = np.real((np.exp(-1j*(PhiOut)*MagOrder) * (c_sens_abs[MagOrder-1])))
    print(B_Main_Rotated)
    
    #Rotation of all components
    c_sens_abs_rot=np.zeros(len(c_sens_abs)).astype('complex')
    c_sens_cmp_rot=np.zeros(len(c_sens_cmp)).astype('complex')
       
    for n in np.arange(len(c_sens_abs)):
        
        c_sens_abs_rot[n]=c_sens_abs[n]*(np.exp(-1j*PhiOut*(n+1)))
        c_sens_cmp_rot[n]=c_sens_cmp[n]*(np.exp(-1j*PhiOut*(n+1)))
        

    #Apply Bucking ratios------------------------------------------------------------

    Buck_ratio=abs(f_abs/f_cmp)

    #Center  localizations
    
    if MagOrder==1:
        Cn_1 = c_sens_cmp_rot[10]
        Cn_2 = c_sens_cmp_rot[11]
        zR = -(Cn_1/(10*Cn_2))
    else:
        Cn_1 = c_sens_abs_rot[MagOrder]
        Cn_2 = c_sens_abs_rot[(MagOrder-1)]
        zR = -(Cn_1/((MagOrder-1)*Cn_2))
        
    DeltaZ=-zR/Rref
    x=np.real(zR)
    y=np.imag(zR)
    
    #Feed Down------------------------------------------------------------
    
    if 'fed'in AnalysisOptions:
    
        #Absolute C's
        
        c_fd_abs=np.zeros(len(c_sens_abs_rot)).astype('complex')   
        
        
        for n in range(1,(len(c_sens_abs_rot)+1)):
            
            
            for k in range(n,(len(c_sens_abs_rot)+1)):
                
                
                c=(math.factorial(k-1)/(math.factorial(n-1)*math.factorial(k-n)))*c_sens_abs_rot[k-1]*(DeltaZ/Rref)**(k-n)
                
                c_fd_abs[n-1]=c_fd_abs[n-1]+c
        
        #Compensated C's    
        
        c_fd_cmp=np.zeros(len(c_sens_cmp_rot)).astype('complex')    
        
        for n in range(1,(len(c_sens_cmp_rot)+1)):
            
            for k in range(n,(len(c_sens_cmp_rot)+1)):
                
                c=(math.factorial(k-1)/(math.factorial(n-1)*math.factorial(k-n)))*c_sens_cmp_rot[k-1]*(DeltaZ/Rref)**(k-n)
                
                c_fd_cmp[n-1]=c_fd_cmp[n-1]+c
    
    else:
        c_fd_abs=c_sens_abs_rot
        c_fd_cmp=c_sens_cmp_rot
        
    #Normalization------------------------------------------------------------

    
    if 'nor'in AnalysisOptions: 
        c_fd_nor_cmp=np.zeros(len(c_fd_cmp)).astype("complex")
        c_fd_nor_abs=np.zeros(len(c_fd_abs)).astype("complex")
        
        
        
        if MagOrder==1:
            #c_fd_nor_cmp[0]=c_fd_cmp[0]  
            c_fd_nor_cmp[0]=B_Main_Rotated
            for m in range(1,len(c_fd_nor_cmp)):     
                c_fd_nor_cmp[m]=10000*(c_fd_cmp[m]/B_Main_Rotated)
        else:
            for n in np.arange(len(c_fd_nor_abs)):
                c_fd_nor_cmp[n]=10000*(c_fd_cmp[n]/B_Main_Rotated)
            
        
        if MagOrder==1:
            c_fd_nor_abs[0]=c_fd_abs[0]                  
            for m in range(1,len(c_fd_nor_abs)):               
                c_fd_nor_abs[m]=10000*(c_fd_abs[m]/B_Main_Rotated)               
        else:
            for n in np.arange(len(c_fd_nor_abs)):
                c_fd_nor_abs[n]=10000*(c_fd_abs[n]/B_Main_Rotated)
    
    
    
    
    
        
    # c_fd_nor_cmp=np.zeros(len(c_fd_cmp)).astype("complex")
  
    # c_fd_nor_cmp[0]=B_Main_Rotated
    # for m in range(1,len(c_fd_nor_cmp)):     
    #     c_fd_nor_cmp[m]=c_fd_cmp[m]
    
    
    #Apply bucking ratio------------------------------------------------------------
    #???
    
    for n in range(MagOrder):
        c_fd_nor_cmp[n]=np.complex(Buck_ratio[n],Buck_ratio[n])
        
    Norm_abs=np.real(c_fd_nor_abs)
    Skew_abs=np.imag(c_fd_nor_abs)
    Norm_cmp=np.real(c_fd_nor_cmp)
    Skew_cmp=np.imag(c_fd_nor_cmp)
    

    
    return [Norm_abs,
            Skew_abs,
            Norm_cmp,
            Skew_cmp,
            x,
            y,
            PhiOut]


def Parameters(path):
    with open(path,'r') as params_file:
        AnalysisOptions=''
        for line in params_file:
            parameter=line.split(':')[0].split('.')[-1]
            if parameter=='samples':
                p_turn=int(line.split(':')[1])
            elif parameter== 'refRadius':
                Rref=float(line.split(':')[1])
            #------------Analysis options---------------
            elif parameter== 'cel':
                AnalysisOptions=AnalysisOptions+' '+parameter
            elif parameter== 'deb':
                AnalysisOptions=AnalysisOptions+' '+parameter
            elif parameter== 'dri':
                AnalysisOptions=AnalysisOptions+' '+parameter
            elif parameter== 'nor':
                AnalysisOptions=AnalysisOptions+' '+parameter
            elif parameter== 'rot':
                AnalysisOptions=AnalysisOptions+' '+parameter

    #p_turn=512
    num_FDI=2
    MagOrder=1
    #Rref=0.050

    return(p_turn,num_FDI,MagOrder,Rref,AnalysisOptions)


def ContinuousRotatingCoilAnalysis (datafile, p_turn, knAbs, knCmp, MagOrder, Rref, AnalysisOptions):
    
    
    df = pd.read_csv(datafile, sep='\s+',names=['Time','df_abs','df_cmp','curr'])
    
    
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
    
    # Multip_list=[NAbs,SAbs,NCmp,SCmp] # Plot multiples Obtained at each turn
    # for dfr in Multip_list:
    #     plt.figure()
    #     for column in dfr.columns[1:]:
    #         plt.plot(dfr["n"],dfr[column])
            
    # Save the dataframes
    
    NAbs.to_excel(datafile.replace("raw_measurement_data","C_Normal_Absolute").replace('.txt','.xlsx'))
    SAbs.to_excel(datafile.replace("raw_measurement_data","C_Skew_Absolute").replace('.txt','.xlsx'))
    NCmp.to_excel(datafile.replace("raw_measurement_data","C_Normal_Compensated").replace('.txt','.xlsx'))
    SCmp.to_excel(datafile.replace("raw_measurement_data","C_Skew_Compensated").replace('.txt','.xlsx'))
    Position.to_excel(datafile.replace("raw_measurement_data","Position").replace('.txt','.xlsx'))
    
    NAbs_ave=pd.DataFrame()
    NAbs_ave['n']=NAbs['n']
    NAbs_ave['NAbs Average']=NAbs[NAbs.columns[1:]].mean(axis=1)
    
    SAbs_ave=pd.DataFrame()
    SAbs_ave['n']=SAbs['n']
    SAbs_ave['SAbs Average']=SAbs[SAbs.columns[1:]].mean(axis=1)
    
    NCmp_ave=pd.DataFrame()
    NCmp_ave['n']=NCmp['n']
    NCmp_ave['NCmp Average']=NCmp[NCmp.columns[1:]].mean(axis=1)
    
    SCmp_ave=pd.DataFrame()
    SCmp_ave['n']=SCmp['n']
    SCmp_ave['SCmp Average']=SCmp[SCmp.columns[1:]].mean(axis=1)
    
    #returns the average of each multipole
    

    return(NAbs_ave,SAbs_ave,NCmp_ave,SCmp_ave)
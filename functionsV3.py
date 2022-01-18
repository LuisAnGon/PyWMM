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
    
    # print(c_sens_abs[0])
    
    
    
    SignalPhase=np.angle(c_sens_abs[(MagOrder-1)])
    # print('Modulus:',np.absolute(c_sens_abs[(MagOrder-1)]))
    # print('Phase:',SignalPhase/np.pi)
    
    if (SignalPhase>np.pi/2):
        # print('SignalPhase>np.pi/2')
        SignalPhase = (SignalPhase-np.pi)
    elif (SignalPhase<-np.pi/2):
        # print('SignalPhase<-np.pi/2')
        SignalPhase = (SignalPhase+np.pi)
    
    
    #direction of the main component
    
    PhiOut=SignalPhase/MagOrder
    # PhiOut=0.437205*np.pi
    #print("PhiOut",PhiOut)
    
    # PhiOut=0.937205*np.pi
    
    
    B_Main_Rotated_norm = np.real((np.exp(-1j*(PhiOut)*MagOrder) * (c_sens_abs[MagOrder-1])))
    B_Main_Rotated_skew = np.imag((np.exp(-1j*(PhiOut)*MagOrder) * (c_sens_abs[MagOrder-1])))
    # print('B_Main_Rotated_norm:',B_Main_Rotated_norm)
    # print('B_Main_Rotated_skew:',B_Main_Rotated_skew)
    
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
        
    # DeltaZ=-zR/Rref
    z=Rref*zR
    x=np.real(z)
    y=np.imag(z)
    
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

    #### FOR THE FILE THAT WE WILL GENERATE WITH PYTHON WE WANT THE NON NORMALIZED VALUES
    #### HENCE WE COMMENT OUT THE ORIGINAL CODE
    
    if 'nor'in AnalysisOptions: 
        c_fd_nor_cmp=np.zeros(len(c_fd_cmp)).astype("complex")
        # c_fd_nor_abs=np.zeros(len(c_fd_abs)).astype("complex")
        
        
        
        if MagOrder==1:
            #c_fd_nor_cmp[0]=c_fd_cmp[0]  
            c_fd_nor_cmp[0]=B_Main_Rotated_norm
            for m in range(1,len(c_fd_nor_cmp)):     
                c_fd_nor_cmp[m]=10000*(c_fd_cmp[m]/B_Main_Rotated_norm)
        else:
            for n in np.arange(len(c_fd_nor_cmp)):
                c_fd_nor_cmp[n]=10000*(c_fd_cmp[n]/B_Main_Rotated_norm)
            
        
        # if MagOrder==1:
        #     c_fd_nor_abs[0]=c_fd_abs[0]                  
        #     for m in range(1,len(c_fd_nor_abs)):               
        #         c_fd_nor_abs[m]=10000*(c_fd_abs[m]/B_Main_Rotated_norm)               
        # else:
        #     for n in np.arange(len(c_fd_nor_abs)):
        #         c_fd_nor_abs[n]=10000*(c_fd_abs[n]/B_Main_Rotated_norm)
    
    
    
    
    
        # This dataframe is for comparing with Roxie
        c_fd_nor_cmp_roxie=np.zeros(len(c_fd_cmp)).astype("complex")
      
        c_fd_nor_cmp_roxie[0]=np.complex(B_Main_Rotated_norm,B_Main_Rotated_skew)
        for m in range(1,len(c_fd_nor_cmp_roxie)):     
            c_fd_nor_cmp_roxie[m]=c_fd_cmp[m]
    
    
    #Apply bucking ratio------------------------------------------------------------
    #???
    
    # for n in range(MagOrder):
    #     c_fd_nor_cmp[n]=np.complex(Buck_ratio[n],Buck_ratio[n])
        
    # Norm_abs=np.real(c_fd_nor_abs)
    # Skew_abs=np.imag(c_fd_nor_abs)
    Norm_cmp=np.real(c_fd_nor_cmp)
    Skew_cmp=np.imag(c_fd_nor_cmp)
    
    
    return [Norm_cmp,
            Skew_cmp,
            x,
            y,
            PhiOut]
    
    # return [Norm_abs,
    #         Skew_abs,
    #         Norm_cmp,
    #         Skew_cmp,
    #         x,
    #         y,
    #         PhiOut]


def Parameters(path):
    
    with open(path,'r') as params_file:
        AnalysisOptions=''
        MagOrder=0
        for line in params_file:
            parameter=line.split(':')[0].split('.')[-1]
            if parameter=='samples':
                p_turn=int(line.split(':')[1])
            elif parameter== 'refRadius':
                Rref=float(line.split(':')[1])
            # elif parameter=="magnetType":
            #     order=str(line.split(':')[1])
            #     if order=="Dipole":
            #         MagOrder=1
            #     if order=="Quadrupole":
            #         MagOrder=2
            
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


def ContinuousRotatingCoilAnalysis (datafile, p_turn, knAbs, knCmp, MagOrder, Rref, AnalysisOptions,step):
    
    
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
    # Position=pd.DataFrame(columns=("Turn","Rref","Options","x","y","Angle"))
    Header=pd.DataFrame()
    
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
        [Norm_cmp,
         Skew_cmp,
         x,
         y,
         PhiOut]=RotatingCoilAnalysisTurn(Fabs, Fcmp, knAbs, knCmp, MagOrder, Rref, AnalysisOptions)
        
        
        
        
        
        # Add the Multipoles corresponding to each turn to different Dataframes depending on their NOrmal/Skew and Absolute/Compensated
        head=[0,Rref,x,y,PhiOut]
        # NAbs["Turn_"+str(i)]=Norm_abs.tolist()
        # SAbs["Turn_"+str(i)]=Skew_abs.tolist()
        
        Header["Turn_"+str(i)]=head
        NCmp["Turn_"+str(i)]=Norm_cmp.tolist()
        SCmp["Turn_"+str(i)]=Skew_cmp.tolist()
        # Position["Turn_"+str(i)]=
        
        
        
        i=i+1
    
    Multip_list=[Header,NCmp,SCmp] # Add number of the Multipole
    for dfr in Multip_list:
        
        dfr=dfr.insert(0,"n",dfr.index+1)
    

    
    # Multip_list=[NAbs,SAbs,NCmp,SCmp] # Plot multiples Obtained at each turn
    # for dfr in Multip_list:
    #     plt.figure()
    #     for column in dfr.columns[1:]:
    #         plt.plot(dfr["n"],dfr[column])
            
    # Save the dataframes
    
    # NAbs.to_excel(datafile.replace("raw_measurement_data","C_Normal_Absolute").replace('.txt','.xlsx'))
    # SAbs.to_excel(datafile.replace("raw_measurement_data","C_Skew_Absolute").replace('.txt','.xlsx'))
    NCmp.to_excel(datafile.replace("raw_measurement_data","C_Normal_Compensated").replace('.txt','.xlsx'))
    SCmp.to_excel(datafile.replace("raw_measurement_data","C_Skew_Compensated").replace('.txt','.xlsx'))
    Header.to_excel(datafile.replace("raw_measurement_data","Header").replace('.txt','.xlsx'))
    
    # NAbs_ave=pd.DataFrame()
    # NAbs_ave['n']=NAbs['n']
    # NAbs_ave['NAbs Average']=NAbs[NAbs.columns[1:]].mean(axis=1)
    
    # SAbs_ave=pd.DataFrame()
    # SAbs_ave['n']=SAbs['n']
    # SAbs_ave['SAbs Average']=SAbs[SAbs.columns[1:]].mean(axis=1)
    
    Header_ave=pd.DataFrame()
    
    Header_ave['Average']=Header[Header.columns[1:]].mean(axis=1)
    Header_ave.index=["position","Rref","x","y","Angle"]
    Header_ave=Header_ave.append(pd.Series(name='Space'))
    
    NCmpindex=["B1"]
    NCmp_ave=pd.DataFrame()
    
    NCmp_ave['Average']=NCmp[NCmp.columns[1:]].mean(axis=1)
    for i in range(len(NCmp)+1)[2:]:
        NCmpindex.append("b"+str(i))
    NCmp_ave.index=NCmpindex
    NCmp_ave=NCmp_ave.append(pd.Series(name='Space'))
    
    
    SCmpindex=["A1"]
    SCmp_ave=pd.DataFrame()
    
    SCmp_ave['Average']=SCmp[SCmp.columns[1:]].mean(axis=1)
    for i in range(len(SCmp)+1)[2:]:
        SCmpindex.append("a"+str(i))
    SCmp_ave.index=SCmpindex
    SCmp_ave=SCmp_ave.append(pd.Series(name='Space'))
    
    
    turn=Header_ave.append(NCmp_ave.append(SCmp_ave))
    
    
    
    
    #returns the average of each multipole
    return turn
    # return(Header,NCmp_ave,SCmp_ave)
    # return(NAbs_ave,SAbs_ave,NCmp_ave,SCmp_ave)
    
    
def ReadFromC(path,step):
    summary=pd.read_csv(path,sep='\t')  
    for i,row in summary.iterrows():
        if (i % 2) != 0:
            summary.drop(i,axis=0,inplace=True)
    
    new_cols=[]
    for col in summary.columns:
        if col != 'Options':
            summary[col]=summary[col].astype('float')
        new_cols.append(col.split('(m)')[0].split('(A)')[0].split(' ')[0])
    
    summary.columns=new_cols
    summary=summary.drop(columns='Unnamed:')
    
       
            
    summarypos=summary[summary['I']>0].reset_index()
    summaryneg=summary[summary['I']<0].reset_index()
    
    
    summarypos['index']=np.array(summarypos.index)*step
    
    summarypos.rename(columns={'index':'position'},inplace=True)
    summarypos=summarypos.transpose()
    summarypos['indicator']=summarypos[summarypos.columns[0]]
    

    for i,row in summarypos.iterrows():
        
        if i[0]=='b':
            summarypos['indicator'][i]='b'
        elif i[0]=='a':
            summarypos['indicator'][i]='a'
        else:
            summarypos['indicator'][i]='h'
                
        bes= summarypos[summarypos['indicator']=='b']      
        aes= summarypos[summarypos['indicator']=='a'] 
        head=summarypos[summarypos['indicator']=='h']

    bes=bes.append(pd.Series(name='Space'))
    aes=aes.append(pd.Series(name='Space'))
    head=head.append(pd.Series(name='Space'))    
    summarypos=head.append(bes).append(aes).drop("indicator",axis=1)
    
    
    summaryneg['index']=np.array(summaryneg.index)*step
    summaryneg.rename(columns={'index':'position'},inplace=True)
    summaryneg=summaryneg.transpose()
    summaryneg['indicator']=summaryneg[summaryneg.columns[0]]
    
    
    for i,row in summaryneg.iterrows():
        
        if i[0]=='b':
            summaryneg['indicator'][i]='b'
        elif i[0]=='a':
            summaryneg['indicator'][i]='a'
        else:
            summaryneg['indicator'][i]='h'
                
        bes= summaryneg[summaryneg['indicator']=='b']      
        aes= summaryneg[summaryneg['indicator']=='a'] 
        head=summaryneg[summaryneg['indicator']=='h']

    bes=bes.append(pd.Series(name='Space'))
    aes=aes.append(pd.Series(name='Space'))
    head=head.append(pd.Series(name='Space'))    
    summaryneg=head.append(bes).append(aes).drop("indicator",axis=1)
    
    summarypos.to_excel('\\'.join(path.split('\\')[0:-1])+'\\summary_pos.xlsx')
    summaryneg.to_excel('\\'.join(path.split('\\')[0:-1])+'\\summary_neg.xlsx')
    
    return(summarypos,summaryneg)


def GetSensitivities(senspath):
    sensitivities = pd.read_csv(senspath, sep='\s+',names=['abs_real','abs_im','cmp_real','cmp_im'])
    knAbs=np.array(sensitivities['abs_real']+ 1j * sensitivities['abs_im'])
    knCmp=np.array(sensitivities['cmp_real']+ 1j * sensitivities['cmp_im'])
    return(knAbs,knCmp)

def Differences(Av,Summ):
    Av_ind=list(Av.index)
    Summ_ind=list(Summ.index)
    coincide1=Summ.index.isin(Av_ind)
    Summ=Summ[coincide1]
    coincide2=Av.index.isin(Summ_ind)
    Av=Av[coincide2]
    
    Av=Av.drop("Space",axis=0).astype("int")
    Summ=Summ.drop("Space",axis=0).astype("int")
    Av=Av.reindex(list(Summ.index))
    # print (Av.index==Summ.index)
    
    
    difference=Av.equals(Summ)

    return difference
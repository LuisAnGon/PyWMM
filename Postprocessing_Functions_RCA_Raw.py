# -*- coding: utf-8 -*-
"""
Created on Tue May 24 12:58:44 2022

@author: luis gonzalez
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Aug  5 13:57:29 2021

THIS FILE CONTAINS

RotatingCoilAnalysisTurn
ContinuousRotatingCoilAnalysis
ReadRoxie
interpolate_Roxie_MM
GetSensitivities
ReadFromCplus
Parameters


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
import scipy
import io
import matplotlib.pyplot as plt
import os
import pyperclip


def ContinuousRotatingCoilAnalysisRaw (df, knAbs, knCmp, MagOrder, Rref, AnalysisOptions,step,skew,rot,CalibAng):
        
    # # =============================================================================
    # #    Takes the fluxes data and puts it in a df
    # # =============================================================================
    #     # df = pd.read_csv(datafile, sep='\s+',names=['Time','df_abs','df_cmp','curr'])
    #     df = pd.read_csv(datafile, sep='\s+',names=['df_abs','time1','df_cmp','time2'])
    # # =============================================================================
    # #    Takes the fluxes data and puts it in a df
    # # =============================================================================
        
        
    # # =============================================================================
    # #     Separates in turns
    # # =============================================================================
    #     turns=len(df)/p_turn
    
    #     turnlst=[]
    #     i=0
    
    #     for turn in np.arange(int(turns)):
        
    #         turn_data=df.iloc[i*p_turn:(i+1)*p_turn].reset_index()
    #         turnlst.append(turn_data)
    #         i=i+1
    # # =============================================================================
    # #     Separates in turns
    # # =============================================================================
        
        
        
        NCmp=pd.DataFrame()
        SCmp=pd.DataFrame()
        Header=pd.DataFrame()  

        df=df.transpose()

        i=0 
        for turn in df.columns:
            # print("*********")
            turn_abs=df[turn].iloc[17:529]
            # print(len(turn_abs))
            turn_abs=turn_abs.reindex(index=np.roll(turn_abs.index,1))#correct data (shift one step - until introduced in FFMM)
            
            turn_cmp=df[turn].iloc[529:1041]
            # print(len(turn_cmp))
            turn_cmp=turn_cmp.reindex(index=np.roll(turn_cmp.index,1))#correct data (shift one step - until introduced in FFMM)
            
            turn_dt=df[turn].iloc[1041:len(df[turn])]
            # print(len(turn_dt))
            turn_dt=turn_dt.reindex(index=np.roll(turn_dt.index,1))#correct data (shift one step - until introduced in FFMM)
            
            
            
            
            if rot<0:
                turn_abs=-turn_abs.iloc[::-1]
                turn_cmp=-turn_cmp.iloc[::-1]
                # turn_abs=turn_abs
                # turn_cmp=turn_cmp
            
    # =============================================================================
    #         Drift Correction of each turn
    # =============================================================================
            if 'dri'in AnalysisOptions:
                
                #Drift Piotr
                
                time=turn_dt.sum()
                sumabs=turn_abs.sum()
                sumcmp=turn_cmp.sum() 
                
                
                Fabs=np.array(turn_abs)-(sumabs/time)*np.array(turn_dt)   
                Fabs=Fabs.cumsum()
                Fcmp=np.array(turn_cmp)-(sumcmp/time)*np.array(turn_dt)
                Fcmp=Fcmp.cumsum()
                
                # Drift Luis
                
                # Fabs = (turn_abs-turn_abs.mean()).cumsum()-turn_abs.cumsum().mean()
                # Fcmp = (turn_cmp-turn_cmp.mean()).cumsum()-turn_cmp.cumsum().mean()
            else:
                
                Fabs=turn_abs.cumsum()
                Fcmp=turn_abs.cumsum()     
    # =============================================================================
    #         Drift Correction of each turn
    # =============================================================================
            
            
    # =============================================================================
    #         Analysis of each turn
    # =============================================================================       
            [Norm_cmp,
             Skew_cmp,
             x,
             y,
             Ang]=RotatingCoilAnalysisTurn(Fabs, Fcmp, knAbs, knCmp, MagOrder, Rref, AnalysisOptions,skew,CalibAng)
            print ("Ang= ",Ang)
    # =============================================================================
    #         Analysis of each turn
    # =============================================================================       
        
        
    # =============================================================================
    #         Creates a df with the shape:
    #         
    #                      Average (Turn1 ... Turn10)
    #             Rref
    #             x
    #             y
    #             phiOut
    #             B1
    #             ...
    #             B15
    #             A1
    #             ...
    #             A15
    # =============================================================================
                
            
            head=[0,Rref,x,y,Ang]        
            Header["Turn_"+str(i)]=head
            NCmp["Turn_"+str(i)]=Norm_cmp.tolist()
            SCmp["Turn_"+str(i)]=Skew_cmp.tolist()

            
            # Position["Turn_"+str(i)]=
            
            i=i+1
            
        
        
        Multip_list=[Header,NCmp,SCmp] # Add number of the Multipole
        
        
        for dfr in Multip_list:   
            dfr=dfr.insert(0,"n",dfr.index+1)
        
        Header_ave=pd.DataFrame()
        
        Header_ave['Average']=Header[Header.columns[1:]].mean(axis=1)
        Header_ave.index=["position","Rref [m]","x [m]","y [m]","Angle (Applied Calib) [mrad]"]
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
    # =============================================================================
    #     Creates a df with the shape:
    # =============================================================================
        
        return turn
    
        
        
def RotatingCoilAnalysisTurn(Fabs, Fcmp, knAbs, knCmp, MagOrder, Rref, AnalysisOptions,skew,CalibAng):
        
    # =============================================================================
    #     Fourier transform
    # =============================================================================
        f_abs=2*fft(np.array(Fabs))/len(Fabs)
        # print("Angle_f_abs: %.64f" %(np.angle(f_abs[(MagOrder-1)])))
        f_cmp=2*fft(np.array(Fcmp))/len(Fcmp)
        f_abs=f_abs[1:16]
        f_cmp=f_cmp[1:16]
    # =============================================================================
    #     Fourier transform
    # =============================================================================
     
    
       
        
    # =============================================================================
    #     Apply coil sensitivity
    # =============================================================================
        Rref_vec=np.zeros(len(f_abs))
        
        for n in np.arange(0,len(Rref_vec)):
            Rref_vec[n]=Rref**n
        
        c_sens_abs = Rref_vec*(1/knAbs)*f_abs
        c_sens_cmp = Rref_vec*(1/knCmp)*f_cmp
    # =============================================================================
    #     Apply coil sensitivity
    # =============================================================================
    
    
    
    
    
    
        
    # =============================================================================
    #     Angle Correction (Apply Rotation)
    # ============================================================================= 
        SignalPhase=np.angle(c_sens_abs[(MagOrder-1)])
        # print("signalPhase= ",SignalPhase)
        
        if (SignalPhase>np.pi/2):
            SignalPhase = (SignalPhase-np.pi)
        
        elif (SignalPhase<-np.pi/2):
            SignalPhase = (SignalPhase+np.pi)
        
        
        #direction of the main component
        
        PhiOut=SignalPhase/MagOrder
        
        
        
        
        if  skew=="SKEW":
            PhiOut=PhiOut-(np.pi/2)
        
        Ang=PhiOut*1000
        Ang=Ang-CalibAng
        
        # print ("PhiOut= ",PhiOut)
        B_Main_Rotated = (np.exp(-1j*(PhiOut)*MagOrder) * (c_sens_abs[MagOrder-1]))
        
        B_Main_Rotated_norm = np.real((np.exp(-1j*(PhiOut)*MagOrder) * (c_sens_abs[MagOrder-1])))
        B_Main_Rotated_skew = np.imag((np.exp(-1j*(PhiOut)*MagOrder) * (c_sens_abs[MagOrder-1])))
        
        if skew=="SKEW":
            B_Main_mod=B_Main_Rotated_skew
        else:
            B_Main_mod=B_Main_Rotated_norm
        
        #Rotation of all components
        c_sens_abs_rot=np.zeros(len(c_sens_abs)).astype('complex')
        c_sens_cmp_rot=np.zeros(len(c_sens_cmp)).astype('complex')
           
        for n in np.arange(len(c_sens_abs)):
            
            c_sens_abs_rot[n]=c_sens_abs[n]*(np.exp(-1j*PhiOut*(n+1)))
            c_sens_cmp_rot[n]=c_sens_cmp[n]*(np.exp(-1j*PhiOut*(n+1)))
    # =============================================================================
    #     Angle Correction (Apply Rotation)
    # ============================================================================= 
    
    
            
    
    
    
    # =============================================================================
    #     Apply Bucking ratios
    # =============================================================================
        Buck_ratio=abs(f_abs/f_cmp)
    # =============================================================================
    #     Apply Bucking ratios
    # =============================================================================
    
    
    
    
    
    # =============================================================================
    #    Center  localizations
    # =============================================================================
        
        if MagOrder==1:
            p=[0,0,0,0,0,0]
            cost=[0,0,0,0,0]
            n=10
            for k in range (n,16):
                p[15-k]=(math.factorial(k-1)/(math.factorial(n-1)*math.factorial(k-n)))*c_sens_cmp_rot[k-1]*(1/Rref)**(k-n)
            r=np.roots(p)
            for i in [0,1,2,3,4]:
                tmp4 = [0]*15
                zR = r[i]/Rref
                for n in range (1,15):
                    for k in range (n,15):
                        tmp3 = (math.factorial(k-1)/(math.factorial(k-n)*math.factorial(n-1)))*c_sens_cmp_rot[k-1]*(zR**(k-n));
                        tmp4[n] = tmp4[n] + tmp3
                
                Cprim = tmp4
                
                for l in [4,5,6,7]:
                    cost[i]=cost[i] + abs(Cprim[2*l])/abs(c_sens_cmp_rot[2*l])
                val=min(cost)
                ind=cost.index(val)
                
                z=r[ind]
                
            
            
                
            
            
            
            # Cn_1 = c_sens_cmp_rot[10]
            # Cn_2 = c_sens_cmp_rot[11]
            # zR = -(Cn_1/(10*Cn_2))
        else:
            Cn_1 = c_sens_abs_rot[MagOrder]
            Cn_2 = c_sens_abs_rot[(MagOrder-1)]
            zR = -(Cn_1/((MagOrder-1)*Cn_2))
            

        x=np.real(z)*1000
        y=np.imag(z)*1000
        
        zR=z/Rref
        
        
        
    # =============================================================================
    #    Center  localizations
    # =============================================================================
        
    
    
    
    
    # =============================================================================
    #     Feed Down
    # =============================================================================   
        if 'fed'in AnalysisOptions:
        
            #Absolute C's
            
            c_fd_abs=np.zeros(len(c_sens_abs_rot)).astype('complex')   
            
            
            for n in range(1,(len(c_sens_abs_rot)+1)):
                
                
                for k in range(n,(len(c_sens_abs_rot)+1)):
                    
                    
                    c=(math.factorial(k-1)/(math.factorial(n-1)*math.factorial(k-n)))*c_sens_abs_rot[k-1]*(zR)**(k-n)
                    
                    c_fd_abs[n-1]=c_fd_abs[n-1]+c
            
            #Compensated C's    
            
            c_fd_cmp=np.zeros(len(c_sens_cmp_rot)).astype('complex')    
            
            for n in range(1,(len(c_sens_cmp_rot)+1)):
                
                for k in range(n,(len(c_sens_cmp_rot)+1)):
                    
                    c=(math.factorial(k-1)/(math.factorial(n-1)*math.factorial(k-n)))*c_sens_cmp_rot[k-1]*(zR)**(k-n)
                    
                    c_fd_cmp[n-1]=c_fd_cmp[n-1]+c
        
        else:
            c_fd_abs=c_sens_abs_rot
            c_fd_cmp=c_sens_cmp_rot
    # =============================================================================
    #     Feed Down
    # =============================================================================
        
    
    
    
    
    # =============================================================================
    #     Normalization
    # =============================================================================   
        if 'nor'in AnalysisOptions: 
            c_fd_nor_cmp=np.zeros(len(c_fd_cmp)).astype("complex")
            # c_fd_nor_abs=np.zeros(len(c_fd_abs)).astype("complex")
            
            
            
            if MagOrder==1:
                c_fd_nor_cmp[0]=B_Main_Rotated
                for m in range(1,len(c_fd_nor_cmp)):     
                    c_fd_nor_cmp[m]=10000*(c_fd_cmp[m]/B_Main_mod)
            else:
                for n in np.arange(len(c_fd_nor_cmp)):
                    c_fd_nor_cmp[n]=10000*(c_fd_cmp[n]/B_Main_mod)
            
            # This dataframe is for comparing with Roxie
            c_fd_nor_cmp_roxie=np.zeros(len(c_fd_cmp)).astype("complex")
          
            c_fd_nor_cmp_roxie[0]=np.complex(B_Main_Rotated_norm,B_Main_Rotated_skew)
            for m in range(1,len(c_fd_nor_cmp_roxie)):     
                c_fd_nor_cmp_roxie[m]=c_fd_cmp[m]
    # =============================================================================
    #     Normalization
    # =============================================================================     
        
    
        Norm_cmp=np.real(c_fd_nor_cmp)
        Skew_cmp=np.imag(c_fd_nor_cmp)
        
        return [Norm_cmp,
                Skew_cmp,
                x,
                y,
                Ang]
    
    
    
def ReadRoxie(file,n=15, NS="NS",bothDipoles=False,norm=False,skew=0):
        """
        Reads Roxie .output in the format that Luis Gonzalez makes the calculations i.e. 30-31 separate plots with B1, B2, B3...A13, A14, A15, (A16-field along path)
    
        Parameters
        ----------
        file : STRING
            DESCRIPTION. Folder where the output file is.
        n" : INT, optional
            DESCRIPTION. number od multipoles The default is 15.
        
        NS : STRING
            DESCRIPTION. Says if the Roxie file contents:
                Only Skew NS=S
                Only Normal NS=N
                Both Normal and Skew NS=NS
                
        "bothDipoles : BOOL, optional
            DESCRIPTION. TRUE if the ROxie results are powering both dipoles, The default is "False.
        "skew : INT, optional
            DESCRIPTION. Says if the magnet geometry is Normal or Skew The default is 0".
    
        Returns
        -------
        list
            [df with the multipoles as a funtion of the longitudinal position,
            df with the  ultipoles at the center]
        """
    # =============================================================================
    #     Takes the Roxie output corresponding to the measured magnet that is in the measurement folder (The Roxie output could be taken from a common folder and chosen according to the specifications given in the measurement folder name)
    # =============================================================================
        for item in os.listdir(file):
            if item.split(".")[-1]=="output":
                fileroxie=item
        f=open(file+"\\"+fileroxie)
        
    # =============================================================================
    #     Takes the Roxie output corresponding to the measured magnet that is in the measurement folder (The Roxie output could be taken from a common folder and chosen according to the specifications given in the measurement folder name)
    # =============================================================================
       
    
        
    
    # =============================================================================
    #     From the .output selects the multipole profile data columns
    # =============================================================================
        s=f.read()
        f.close()
        s=s.replace('/','')
        s=s.replace('NUMBER OF OBJECTIVES AND CONSTRAINTS','GRAPH')
        a=s.split('GRAPH')
        dflst=[]
    # =============================================================================
    #     From the .output selects the multipole profile data columns
    # =============================================================================
    
    
    
    
     
    # =============================================================================
    #     Creates a df with all multipole profiles (alles)
    # =============================================================================
        alles=pd.DataFrame()
        malles=pd.DataFrame()
        c=0
        for i in range(2,(len(a)-1),1):
            
            b=a[i].split('\n')[2:-3]
            
            df=pd.DataFrame(b,columns=['a'])
            dff=df.a.str.split('   ',expand=True)
            dff.columns=['a','b','c','d','e']
            new=pd.DataFrame()
            
            new['position']=dff['d'].astype(float)
            new['Signal']=dff['e'].astype(float)   
            
            
            if NS == "NS":
            
                if c<n:
                    
                    alles['position']=new['position']
                    alles['B'+str(c+1)+" Roxie"]=new['Signal'].fillna(method='ffill')
                    
                else:
                    
                    alles['A'+str(c-(n-1))+" Roxie"]=new['Signal'].fillna(method='ffill')
                    
            elif NS == "N":
                alles['position']=new['position']
                alles['B'+str(c+1)+" Roxie"]=new['Signal'].fillna(method='ffill')
                
            elif NS == "S":
                alles['position']=new['position']
                alles['A'+str(c+1)+" Roxie"]=new['Signal'].fillna(method='ffill')
                # malles['B'+str(odds[c])]=new['Signal']
                # alles['B'+str(odds[c])+'Shift_up']=new['Signal'].shift(periods=1,axis=0,fill_value=new['Signal'][0])  
                # alles['B'+str(odds[c])+'Shift_down']=new['Signal'].shift(periods=-1,axis=0,fill_value=new['Signal'][0])  
            c+=1
    # =============================================================================
    #     Creates a df with all multipole profiles (alles)
    # =============================================================================
    
        
    
       
    # =============================================================================
    # Shifts the position columns for the profile to be centered in 0
    # =============================================================================
        shift=abs(abs(alles['position'].min())-abs(alles["position"].max()))/2
        alles["position"]=alles['position']+shift
        
        alles.to_excel(file+"\\_Alles.xlsx")
    # =============================================================================
    # Shifts the position columns for the profile to be centered in 0
    # =============================================================================
        
    
    
    
    
    # =============================================================================
    # Creates a df with the Multipoles in the center of the magnet (MP) already normalized
    # =============================================================================
        if bothDipoles==False:
    
            MP=10000*(alles.iloc[int(max(alles.index)/2)]/(max(alles.iloc[int(max(alles.index)/2)])))
            MP=MP.drop(['position'],axis=0)        
    # =============================================================================
    # Creates a df with the Multipoles in the center of the magnet (MP) already normalized
    # =============================================================================
    
    
    
    
        correction=10000/max(alles.iloc[int(max(alles.index)/2)])
        colstocorrect=list(alles.columns)
        colstocorrect.remove("position")
        alles[colstocorrect]=alles[colstocorrect]*correction
        
        
    
    
    
        
    
    
    # =============================================================================
    # Normalizes all the profiles (except the main field) to the main field A1 or B1 depending on whether the magnet is Normal or Skew. Then normalizes the Main Field to its maximum
    # =============================================================================
        
        if norm:
            print("NORM")
            if skew==1:
                mainRoxie="A1 Roxie"
            elif skew==0:
                mainRoxie="B1 Roxie"
                
            colstonorm=list(alles.columns)
            colstonorm.remove("position")
            colstonorm.remove(mainRoxie)
            
            for col in colstonorm:
                  alles[col]=10000*alles[col]/alles[mainRoxie]
           
            alles[mainRoxie]=10000*alles[mainRoxie]/(alles[mainRoxie][len(alles[mainRoxie])/2])
    # =============================================================================
    # Normalizes all the profiles (except the main field) to the main field A1 or B1 depending on whether the magnet is Normal or Skew. Then normalizes the Main Field to its maximum
    # =============================================================================
    
    
    
    
    
    # =============================================================================
    # Saves both df´s            
    # =============================================================================
        filename=fileroxie.split(".")[0]
        alles.to_excel(file+"\\"+filename+"_Roxie_MP.xlsx")
        MP.to_excel(file+"\\"+filename+"_Roxie_MP_Centre.xlsx")
    # =============================================================================
    # Saves both df´s 
    # =============================================================================
        
        return [MP,alles]
    
def SelectRoxie(iron,inner,skew):
        """
        Creates the path to input into ReadRoxie according to the variables iron, inner skew
    
        Parameters
        ----------
        iron : int
            0 or 1.
        inner : int
            0 or 1.
        skew : int
            0 or 1.
    
        Returns
        -------
        path : STR
            The path of the corresponding Roxie File
        """
    # =============================================================================
    #     Creates the string part of the file corresponding to the iron - Inner/Outer - orientation Respectively
    # =============================================================================
        if iron==1:
            ironstr="Iron"
        else:
            ironstr="No Iron"
    
        if inner==1:
            innerstr="Inner Dipole Short"
        else:
            innerstr="Outer Dipole"
    
        if skew==1:
            skewstr="Rot90"
        else:
            skewstr="Norm"
    # =============================================================================
    #     Assembles the path string
    # =============================================================================
        path=r"C:\PyWMM\ROXIE-MCBXFB"+"\\"+ironstr+"\\"+ironstr+" "+innerstr+r"\MCBXFB"+" "+ironstr+" "+innerstr+" "+skewstr
    # =============================================================================
        return path
    
def SelectRoxieManual():
        """
        Creates the path to input into ReadRoxie according to the variables iron, inner skew that are inserted MANUALLY by the user
    
        Returns
        -------
        path : str
            The path of the corresponding Roxie File
    
        """
        
        
        
    # =============================================================================
    #     Creates the string part of the file corresponding to the iron - Inner/Outer - orientation Respectively
    # =============================================================================
        iron=input("Iron? (Y/N)").upper()
        if iron=="Y":
            ironstr="Iron"
        elif iron=="N":
            ironstr="NO Iron"
        else:
            print("Answer must be Y or N")
            
        inner=input("Inner or Outer?").upper()
        if inner=="INNER":
            innerstr="Inner Dipole Short"
        elif inner=="OUTER":
            innerstr="Outer Dipole"
        else:
            print("Answer must be Inner or Outer")
            
        skew=input("Normal or Skew?").upper()
        if skew=="NORMAL":
            skewstr="Norm"
        elif skew=="SKEW":
            skewstr="Rot90"
        else:
            print("Answer must be Normal or Skew")
    # =============================================================================
    
    
    # =============================================================================
    #     Assembles the path string
    # =============================================================================
        path=r"C:\ROXIE-MCBXFB"+"\\"+ironstr+"\\"+ironstr+" "+innerstr+r"\MCBXFB"+" "+ironstr+" "+innerstr+" "+skewstr
    # =============================================================================
    
        return path
    
def interpolate_Roxie_MM(folder,dfRoxie,dfMM,name,tipo,save=False):
        """
        Takes the Roxie dataframe in the shape "Alles" and interpolates it to have values in the same positions as the magnetic measurements.
    
        Parameters
        ----------
        folder : string
            route of the folder where you want to save it (in case you want to save it).
        dfRoxie : DataFrame
            It has the shape Alles i.e.
            
            Position | B1 Roxie | ... | B15 Roxie | ... | A1 Roxie | ... | A15 Roxie
            -1500
              ...
             1500
     
        dfMM : Dataframe
            It has the shape of AV_NN or AV_NS.
            
                      |   |   |
            position
            B1
            b2
            ...
            b15
            
                              |   |   |
            position
            A1
            a2
            ...
            a15
            
        name : string
            Name with which the file would be saved (in case you want to save it).
    
        Returns
        -------
        result : DataFrame
            The Roxie calculations interpolated to the same positions of the magnetic measurements.
    
        """
        
    # =============================================================================
    #     Creates a df with one columns which is the positions of the magnetic measurements
    # =============================================================================
        rango_mm=pd.DataFrame(dfMM['position']).astype(float)
        rango_mm.columns=['position']
    
    # =============================================================================
    #     Creates a df with one columns which is the positions of the magnetic measurements
    # =============================================================================
      
    
     
    # =============================================================================
    #     Interpolates Roxie df to have step of 1mm
    # =============================================================================
        dfRoxie=dfRoxie.astype(float)
    
        rango_RX=pd.DataFrame(np.arange(int(dfRoxie['position'].min()),int(dfRoxie['position'].max())+1,1),columns=['position']).astype(float)
    
        result=dfRoxie.merge(rango_RX,how="outer").sort_values(by='position').interpolate(method='linear')
    
    # =============================================================================
    #     Interpolates Roxie df to have step of 1mm
    # =============================================================================
    
    
    
    # =============================================================================
    #     Merge MM and Roxie DFs Eliminating from Roxie df the positions that are not in the MM df
    # =============================================================================   
    
        interpolated=pd.DataFrame(columns=result.columns)
        for pos in rango_mm.position:
            
            section=result[result["position"]>(int(pos)-300)]
            section=section[section["position"]<pos+300]
            row=[section[col].mean() for col in section.columns]
            interpolated.loc[len(interpolated)]=row
        interpolated.position=list(rango_mm.position)
        result=interpolated
        #result.to_excel(file+"//roxie_int_"+name+".xlsx")
        
        
    # =============================================================================
    #     Eliminate from Roxie df the positions that are not in the MM df
    # =============================================================================
    
        return result
    
def GetSensitivities(senspath):
        """
        Reads the coil sensitivities file
    
        Parameters
        ----------
        senspath : String
            path where the sesitivity file is stored
    
        Returns
        -------
        Sensitivities Absolute signal
        Sensitivities Compensated Signal.
    
        """
        sensitivities = pd.read_csv(senspath, sep='\s+',names=['abs_real','abs_im','cmp_real','cmp_im'])
        knAbs=np.array(sensitivities['abs_real']+ 1j * sensitivities['abs_im'])
        knCmp=np.array(sensitivities['cmp_real']+ 1j * sensitivities['cmp_im'])
        return(knAbs,knCmp)
    
def ReadFromCplus(path,step):
        """
        Reads the data generated by the C++ post-processing software
    
        Parameters
        ----------
        path : STRING
            Psth of the C++ file.
        step : FLOAT
            Longitudinal step between measurements
    
        Returns
        -------
        Summary of the results for positive current
        Summary of the results for negative current
    
        """
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
    
def Parameters(path):
        """
        
        Reads the measurement parameters exported by FFMM.
        
        
        Parameters
        ----------
        path : STRING
            Path of the parameters file
    
        Returns
    
        -------
                points per turn
                mag order
                num FDI
                refRadius
                Analysisoptions
    
        """
        
        
        
    # =============================================================================
    #     Iterates line by line looking for the desired values to return
    # =============================================================================
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
                
                #------------Analysis op options---------------
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
        
    # =============================================================================
    #     Iterates line by line looking for the desired values to return
    # =============================================================================
        
        return(p_turn,num_FDI,MagOrder,Rref,AnalysisOptions)
    
    

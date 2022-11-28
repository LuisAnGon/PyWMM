# -*- coding: utf-8 -*-
"""
THIS FILE CONTAINS

Interpolate
shiftit
simulate_coil
Differences

"""

import pandas as pd
import scipy
import io
import matplotlib.pyplot as plt
import numpy as np
import tkinter as tk
from tkinter import filedialog
import os
import pyperclip
import win32clipboard






def interpolate(df,step=0.5):
    """
    This Funtion is meant to interpolate Roxieall-type of dataframes.
    However it can be used with any dataframe.
    The shape shall be:
        
        position |B1|B2|...|B15|A1|...|A15|
        
        ......... Values .................
    
        It always interpolates on the column named "position"

    Parameters
    ----------
    df : DataFrame
    
        position |B1|B2|...|B15|A1|...|A15|
        ......... Values ..................
    
    step : float, optional
        The step of the interpolation. The default is 0.5.

    Returns
    -------
    result : DataFrame
        The interpolated df with the same shape as the original df.

    """
    
    rango=pd.DataFrame(np.arange(int(df['position'].min()),int(df['position'].max())+step,step),columns=['position'])
    rango.position=rango.position.astype("float")
    result=df.merge(rango,how="outer").sort_values(by='position').interpolate(method="linear")
    result=result.merge(rango,how="right")
    
    return result


def shiftit(file,df,shifting):
    """
    Takes a Roxieall df and shifts the spectrum in the longitudinal position

    Parameters
    ----------
    file : STR
        Name of the .xlsx file that will be saved.
    df : Data-Frame
        It is a Roxieall Dataframe.
    shifting : float
        Length in mm that we want the spectrum to be shifted.

    Returns
    -------
    dfs : DataFrame
        Shifted DataFrame.

    """
# =============================================================================
#     Defines the step of the Roxie File
# =============================================================================
    step=abs(df.position[2]-df.position[1])
    print("Pace/step: ",step)
# =============================================================================


# =============================================================================
#     Defines the period (number of vertical cells up/down each of the df columns) that the data in each column must be shifted
# =============================================================================
    period=int(shifting/step)
    print("Period: ",period)
# =============================================================================



# =============================================================================
# Creates new DataFrame and with new columns anad applies the df.shift() method with the period previously defined    
# =============================================================================
    dfs=pd.DataFrame()
    if period<1:
        print('Error: Shfting (in mm) must be larger than the data step: '+str(period)+' mm')
        return None
    else:
        dfs['position']=df['position']
        for col in df.columns[1:]:
            # print(col)
            dfs[str(col)]= df[col].fillna(method='ffill')
            dfs[str(col)+' shift up']= df[col].shift(periods=period,axis=0,fill_value=df[col][0]).fillna(method='ffill')
            dfs[str(col)+' shift down']= df[col].shift(periods=-period,axis=0,fill_value=df[col][len(df[col])-1]).fillna(method='ffill')
            dfs.to_excel(file+"Multipoles_shifted_"+str(shifting)+"mm.xlsx")
# =============================================================================
    
    return dfs



def simulate_coil(df,coillength,centered=True,displacement=0):
    """
    Simulates the result obtained by a roatating coil
    
    NEEDS - The Multipole profiles FROM ROXIE INTERPOLATED 1mm 
    
    ATENTION - Returns the SUM of all the values corresponding to the RotCoil NOT THE INTEGRAL
    
    It can be chosen to have the coil displaced

    Parameters
    ----------
    df : DataFrame
        It has the shape of Roxieall.
    coillength : INT
        The length of the coil.
    centered : BOOL, optional
        Centered= TRUE means that the center of the magnet is coincident with the center of the coil. 
        Centered=TRUE
        _______________________________________________________________________
        |_________________________________|_____________________________________| <-MAGNET
                                    |___COIL____|                 
        
        Centered= FALSE means thatthe center of the magnet is coincident with the extremity of the coil.
        Centered=FALSE
         _______________________________________________________________________
        |_________________________________|_____________________________________| <-MAGNET
                             |___COIL____|  
        
        The default is True.
    
    displacement : INT, optional
        Used to simulate displacements on the rotating coil. The default is 0.

    Returns
    -------
    integrated : DataFrame
        The simulated coil measurement.

    """
# =============================================================================
#     Define the positions of the coil according to  the coil length and the centering
# =============================================================================   
    integrated=pd.DataFrame(columns=df.columns)
    newpos=[]
    if centered:
        pos=0
        while pos < df.position.max()-coillength:
            pos=pos+coillength
            newpos.append(pos)
            newpos.append(-pos)
        newpos.append(0)
        newpos.sort()
# =============================================================================
#     Define the positions of the coil according to  the coil length and the centering
# =============================================================================



# =============================================================================
#     Apply displacement
# =============================================================================
    newpos=list(map(lambda x: x+displacement, newpos))
# =============================================================================
#     Apply displacement
# =============================================================================


# =============================================================================
#     Sum the values corresponding to each Rotating coil section
# =============================================================================
    df.index=df.position
    df=df.rename_axis("index")
    i=0
    for item in newpos:
        start=item-coillength/2
        end=item+coillength/2
        # print(start, "--", item, "--", end)
        # print(df.loc[start:end].sum())
        integrated=integrated.append(df.loc[start:end].sum()/coillength,ignore_index=True)
        integrated.loc[i,"position"]=item
        # print("*********",i,"**********")
        i=i+1
# =============================================================================
#     Sum the values corresponding to each Rotating coil section
# =============================================================================   
    return integrated


def ReadClipboard(multip=[],norm=False):
    """
    Creates a DataFrame with the copied data from Excel. The data is certain number of multipole profiles copied without column names
    The copied data must have the shape (n columns):
        
        |Position|Multipole(B1,B2...B15...A15)|
        | Values |            Values           |
        

    Parameters
    ----------
    multip : list
        list of strings with the multipoles copied - these will be the column names.
    norm : BOOL, optional
        If TRUE, normalizes a column selected by the user through the console to its value in the middle point (x=0). The default is False.

    Returns
    -------
    df : DataFrame
        A DataFreame with the copied data.

    """
    
# =============================================================================
#   asign a variable to the copied string
# =============================================================================
    a=pyperclip.paste()
# =============================================================================
    


# =============================================================================
#   Split the string into rows
# =============================================================================
    a=a.split("\r\n")[:-1] 
# =============================================================================



# =============================================================================
#   Creates a list of rows, each rows is a list with the values in that row 
# =============================================================================
    rows=[]
    for item in a:
        row=list(item.split("\t"))
        rows.append(row)
# =============================================================================


# =============================================================================
#   Put the rows into a dataframe        
# =============================================================================
    df=pd.DataFrame(data=rows)
# =============================================================================
    
    
    
# =============================================================================
#   Add column names, i.e. the multipoles indicated when calling the function OR generic col-n if no multipoles were indicated
# =============================================================================
    if multip !=[]:
        multip=["position"]+multip
    
    else:
        
        i=0
        for m in range(len(list(df.columns))):
            multip.append("col-"+str(i))
            i+=1
        print (multip)
        print (list(df.columns))
    df.columns=multip
    df=df.astype(float)
# =============================================================================



# =============================================================================
#     Normalizes, if wanted, a column selected by the user through the console to its value in the middle point (x=0)
# =============================================================================
    if norm:
        print(list(df.columns))
        colnorm=input("Column to normalize to?")
        for col in list(df.columns)[1:]:
            if col==colnorm:
                df[colnorm]=df[colnorm]/df[colnorm][int(max(df.index)/2)]
            else:
                df[col]=df[col]*df[colnorm]
# =============================================================================  
    return df




def Differences(Av,Summ):
    """
    ????????????????????????????????
    ????????????????????????????????
    ????????????????????????????????

    Parameters
    ----------
    Av : TYPE
        DESCRIPTION.
    Summ : TYPE
        DESCRIPTION.

    Returns
    -------
    difference : TYPE
        DESCRIPTION.

    """
    Av_ind=list(Av.index)
    Summ_ind=list(Summ.index)
    coincide1=Summ.index.isin(Av_ind)
    Summ=Summ[coincide1]
    coincide2=Av.index.isin(Summ_ind)
    Av=Av[coincide2]
    
    Av=Av.drop("Space",axis=0).astype("int")
    Summ=Summ.drop("Space",axis=0).astype("int")
    Av=Av.reindex(list(Summ.index))

    
    
    difference=Av.equals(Summ)

    return difference

def ClipboardRoxie():
    """
    When the values of the multipoles have been copied from the output file, 
    run ClipboardRoxie() -> The content of the clipboard will be formatted 
    so that a simple paste into Excel will give you two columns

    Returns
    -------
    None.

    """

# get clipboard data
    win32clipboard.OpenClipboard()
    data = win32clipboard.GetClipboardData().replace("\r\n","").replace(" ","").replace("b","\r\nb").replace("a","\r\na").replace(":","\t")
    win32clipboard.CloseClipboard()
    
    
    win32clipboard.OpenClipboard()
    win32clipboard.EmptyClipboard()
    win32clipboard.SetClipboardText(data)
    win32clipboard.CloseClipboard()


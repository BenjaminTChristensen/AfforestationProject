# -*- coding: utf-8 -*-
"""
Created on Thu Nov 18 14:34:02 2021

@author: Henrik Thomsen
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

Rain = pd.read_csv('DMISamletDataFrame.csv', usecols = ['Tid', 'Precipitation'])
ForestType = 'Grass' #Use 'Pine', 'BLT' (Broad-Leaved Trees), 'Grass'
crop = np.ones((12,1))
LV = np.arange(0.5,6.275, 0.275).reshape(21,1)
Pine = np.linspace(1,8,21).reshape(21,1)
LAI = np.ones((21,12)) #Reference LAI matrix
C = 0.005 #cm interception; * LAI later
if ForestType == 'Pine':
    LAI *= Pine
    crop*=1.4
    crop[3:9] = 1.5

elif ForestType =='BLT':
    LAI *= 0.5
    LAI[:,3:9] = LV
    crop*=0.85
    crop[3:9] = 1.05

elif ForestType == 'Grass':
    LAI = np.ones(366)*1.5
    LAI[100:130] = np.linspace(1.5,5,30)
    LAI[130:151] = 5
    LAI[151:165] = np.linspace(1.5, 5, 14)
    LAI[165:188] = 5
    LAI[188:200] = np.linspace(1.5,5,12)
    LAI[200:223] = 5
    LAI[223:250] = np.linspace(1.5,5,27)
    LAI[250:284] = 5
    crop*=1.05

else:
    print("Wrong value in 'ForestType', use 'Pine', 'BLT' or 'Grass'")
#ET0 = pd.read_csv('penman.csv',index_col = 0) #mm/h
#t = 365 * 12 * 20*24
#Rain = pd.read_csv('Rain_hourly.csv')#µm/s

Rain['Precipitation']*=60/10000 # cm/h


Rain['Interception']= C 
"""
LAI_LV = np.array([[0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5],[0.5,0.5,0.5,0.775,0.775,0.775,0.775,0.775,0.775,0.5,0.5,0.5],
                   [0.5,0.5,0.5,1.05,1.05,1.05,1.05,1.05,1.05,0.5,0.5,0.5],[0.5,0.5,0.5,1.325,1.325,1.325,1.325,1.325,1.325,0.5,0.5,0.5],
                   [0.5,0.5,0.5,1.6,1.6,1.6,1.6,1.6,1.6,0.5,0.5,0.5],[0.5,0.5,0.5,1.875,1.875,1.875,1.875,1.875,1.875,0.5,0.5,0.5],
                   [0.5,0.5,0.5,2.15,2.15,2.15,2.15,2.15,2.15,0.5,0.5,0.5],[0.5,0.5,0.5,2.425,2.425,2.425,2.425,2.425,2.425,0.5,0.5,0.5],
                   [0.5,0.5,0.5,2.7,2.7,2.7,2.7,2.7,2.7,0.5,0.5,0.5],[0.5,0.5,0.5,2.975,2.975,2.975,2.975,2.975,2.975,0.5,0.5,0.5],
                   [0.5,0.5,0.5,3.25,3.25,3.25,3.25,3.25,3.25,0.5,0.5,0.5],[0.5,0.5,0.5,3.525,3.525,3.525,3.525,3.525,3.525,0.5,0.5,0.5],
                   [0.5,0.5,0.5,3.8,3.8,3.8,3.8,3.8,3.8,0.5,0.5,0.5],[0.5,0.5,0.5,4.075,4.075,4.075,4.075,4.075,4.075,0.5,0.5,0.5],
                   [0.5,0.5,0.5,4.35,4.35,4.35,4.35,4.35,4.35,0.5,0.5,0.5],[0.5,0.5,0.5,4.625,4.625,4.625,4.625,4.625,4.625,0.5,0.5,0.5],
                   [0.5,0.5,0.5,4.9,4.9,4.9,4.9,4.9,4.9,0.5,0.5,0.5],[0.5,0.5,0.5,5.175,5.175,5.175,5.175,5.175,5.175,0.5,0.5,0.5],
                   [0.5,0.5,0.5,5.45,5.45,5.45,5.45,5.45,5.45,0.5,0.5,0.5],[0.5,0.5,0.5,5.725,5.725,5.725,5.725,5.725,5.725,0.5,0.5,0.5],
                    [0.5,0.5,0.5,6,6,6,6,6,6,0.5,0.5,0.5]])
"""



Rain['Tid'] = pd.to_datetime(Rain['Tid'],format="%Y-%m-%d %H:%M")
Rain['month']=Rain.Tid.dt.month-1
Rain['day'] = Rain.Tid.dt.dayofyear-1
Rain['year'] = Rain.Tid.dt.year-2017
Interep = LAI*C
#Calculating net Precipitation
for i in range(1, len(Rain)):
    if Rain.iloc[i,1] >0 and Rain.iloc[i-1,1] == 0:
        Rain.iloc[i,2] = C
    else:
        Rain.iloc[i,2] = 0

Net=[]
Bru = []
Inter = []
indexer = []
for period in range(8):
    for i in range(len(Rain)):
        Brutto = Rain['Precipitation'][(period*len(Rain)+i)%len(Rain)]
        if ForestType != 'Grass':
            if (period*4+Rain['year'][i]) >= 20:
                Interception = Rain['Interception'][i]*LAI[-1,Rain['month'][i]]**2
            else:
                Interception = Rain['Interception'][i]*LAI[(period*4+Rain['year'][i]),Rain['month'][i]]**2
        else:
            Interception = Rain['Interception'][i]*LAI[Rain['day'][i]]**2
        Netto = Brutto-Interception
        Net.append(max(Netto, 0))
        Bru.append(Brutto)
        Inter.append(Interception)
        indexer.append([period*len(Rain)+i, period, i])
fig, ax = plt.subplots(3,1)
ax[0].plot(Net)
ax[1].plot(Bru)
ax[2].plot(Inter)         
Inter = pd.DataFrame(Inter)
Rain = pd.concat([Rain]*8, ignore_index = True)
Rain['Interception'] = Inter
Overload = False #Determination of the interception at a given rain event exceeds the rain amount
WB = 0 #Initial waterbalance loss, due to no adjusted interception
for i in range(len(Rain)):
    if Overload == True:
        Rain['Interception'][i] = Remaining
    if Rain['Precipitation'][i] == 0 and Rain['Interception'][i] > 0:
        WB = WB + Rain['Interception'][i] #Adding the removed interception, since no rain can be intercepted.
        Rain['Interception'][i] = 0
    if Rain['Interception'][i] > Rain['Precipitation'][i]:
        Remaining = Rain['Interception'][i]-Rain['Precipitation'][i]
        Rain['Interception'][i] = Rain['Precipitation'][i]
        Overload = True
        
    else:
        Overload = False
    
        

print(Inter.sum()) #Checking that the original Interception equals the newly calculated interception + WaterBalance loss
print(Rain['Interception'].sum()+WB)
print(len(Rain.query('Precipitation == 0 and Interception > 0'))) #Counting non-raining events with interception

Rain['Interception'].to_csv('Interception_'+ForestType+'.csv')

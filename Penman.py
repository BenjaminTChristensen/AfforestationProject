# -*- coding: utf-8 -*-
"""
'Created on' %(date)s

@author: %Benjamin Terndrup Christensen
"""

'''Import af diverse packages, som har egenskaber der kan benyttes.
Numpy er god til arrays (vektorer, matricer etc.), tilsvarende med pi, e osv.
Pandas er god til databehandling og datamanipulation af tabeller, blandt andet 
DataFrames og pivot tabeller, samt at læse csv-filer.
Matplotlib bruges til at plotte grafer/figurer.
'''
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime
begin_time = datetime.now()
ForestType = 'BLT' #Use 'Pine', 'BLT' (Broad-Leaved Trees), 'Grass'
crop = np.ones((12,1))
LV = np.arange(0.5,6.275, 0.275).reshape(21,1)
Pine = np.linspace(1,8,21).reshape(21,1)
LAI = np.ones((21,12)) #Reference LAI matrix

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
    LAI = np.ones(365)*1.5
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

#De web-scrapede data fra DMI, Aalborg hentes i en DataFrame: 
#Humid [%], Rad [W/m²], Temp [°C], Wind [m/s], Pressure [hPa]
DMI = pd.read_csv('DMISamletDataFrame.csv')
DMI['Tid'] = pd.to_datetime(DMI['Tid'], format = "%Y-%m-%d %H:%M")
DMI = DMI[DMI['Tid'] < '2021']
DMI['Wind'] *= 4.87/(np.log(67.8*10-5.42)) #Recalculation of v_wind at 2 m, instead of the measured speed at 10 m.
"""Start på udregning af Penmans, først defineres konstanterne:"""
def gamma(P, Cp = 1.013e-3, epsilon = 0.622, L_heat = 2.45):
    a = Cp/(epsilon*L_heat)*P
    return a

def radMJ(rad): #Omregning af W/m² -> MJ/m²*h
    b = rad*3600*1e-6
    return b

Cn = 66
#Omregn Rad fra W/m² til MJ/(m²*hr)

DMI['Radiation'] = radMJ(DMI['Radiation'])

#Lufttryk fra hPa til kPa
DMI['Pressure'] /= 10

#Dataframe med variablerne der regnes.
varss = pd.DataFrame()
#Delta fra Penman regnes
varss['Delta'] = 2503*(0.6108*np.exp(17.27*DMI['TempAvg']/(DMI['TempAvg']+237.3)))/(DMI['TempAvg']+237.3)**2

#Tjekker om parameteren er i dataframen, hvis dette ikke er tilfældet, antages en konstant.

#def Penman(veg, c, k, d):
varss['Cd'] = (DMI['Radiation'] > 0) * (-1.45)+1.7
varss['gamma'] = gamma(DMI['Pressure'])
varss['DT'] = varss['Delta']/(varss['Delta']+varss['gamma']*(1+varss['Cd']*DMI['Wind']))
varss['PT'] = varss['DT'] * varss['gamma']/varss['Delta']
varss['TT'] = (Cn/(DMI['TempAvg']+273))*DMI['Wind']
varss['e_s'] = 0.6108*np.exp(17.27*DMI['TempAvg']/(DMI['TempAvg']+273.3))
varss['e_a'] = DMI['Humidity']/100*varss['e_s']
varss['KG'] = (DMI['Radiation'] > 0)*(-1.6)+2
varss['month'] = DMI.Tid.dt.month-1
varss['year'] = DMI.Tid.dt.year-2017
varss['day'] = DMI.Tid.dt.day-1
varss['ETwind'] = varss['PT']*varss['TT']*(varss['e_s']-varss['e_a'])
if ForestType == 'Grass':
    for i in range(6):
        if i == 5:
            varss['ETrad6'] = varss['DT'] * 0.408 * (DMI['Radiation'])*DMI['Radiation']
        else:
            varss['ETrad'+str(i+1)] = varss['DT'] * 0.408 * (DMI['Radiation']-(varss['KG']*np.exp(-0.5*LAI[varss['day']])*DMI['Radiation']))
        DMI['ET0_'+str(i*4+2017)+'-'+str(i*4+2020)] = (varss['ETrad'+str(i+1)] + varss['ETwind'])/10
else:
    for i in range(6):
        if i == 5:
            varss['ETrad6'] = varss['DT'] * 0.408 * (DMI['Radiation']-(varss['KG']*np.exp(-0.5*LAI[-1,varss['month']])*DMI['Radiation']))
        else:
            varss['ETrad'+str(i+1)] = varss['DT'] * 0.408 * (DMI['Radiation']-(varss['KG']*np.exp(-0.5*LAI[i*4+varss['year'],varss['month']])*DMI['Radiation']))
        DMI['ET0_'+str(i*4+2017)+'-'+str(i*4+2020)] = (varss['ETrad'+str(i+1)] + varss['ETwind'])/10

DMI = DMI.set_index(pd.DatetimeIndex(DMI['Tid']), inplace = False).iloc[:,1:]
num = DMI._get_numeric_data()
num[num<0] = 0  
DMI.iloc[:,-6:].to_csv('Penman_'+ForestType+'.csv')
print(datetime.now()-begin_time)

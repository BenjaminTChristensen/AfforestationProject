import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from datetime import datetime
begin_time = datetime.now()
print(begin_time)
regn = pd.read_csv('DMISamletDataFrame.csv', usecols=['Tid', 'Precipitation'])
regn['Precipitation'] = regn['Precipitation']*60/10000  # From µm/s to cm/h

ForestType = 'BLT'  # Use 'Pine', 'BLT' (Broad-Leaved Trees), 'Grass'
Interception = pd.read_csv('Interception_'+ForestType+'.csv', index_col=0)
Interception = np.array(Interception).reshape(len(Interception))

"""Penmans loaded from the csv file and increasing for every len of the original dataframe.
First ET column is 2017-2020 (DMI data), the following columns are then 2021-2024, 2025-2028 ... 2037+
resulting in an updated ET0, based on the growth of the forest and their LAI"""

ET0 = pd.read_csv('Penman_'+ForestType +'.csv')  # Initial Evapotranspiration. Penmans used later
ET20yr = ET0.iloc[:, -1]  # Fullgrown forest with a constant LAI
ET0 = pd.concat([ET0.iloc[:,1], ET0.iloc[:,2], ET0.iloc[:,3], ET0.iloc[:,4], ET0.iloc[:,5], ET0.iloc[:,6], ET20yr, ET20yr], ignore_index = True)
MB_check =[]


regn['Tid'] = pd.to_datetime(regn['Tid'], format='%Y-%m-%d %H:%M')
mon = regn.Tid.dt.month-1
yr = regn.Tid.dt.year-2017
regn = np.array(regn['Precipitation'])

PM = []
theta_PMMB = []
"""Sensitivity analysis parameter change"""
ParameterChange = 'OG_values' #Write the belowmentioned parameter that is changed
"""Soil properties"""
theta_s = 0.47  # Water content at saturation [cm^3 H2O / cm^3 soil]
K_s = 2.58  # 2.58 Saturated hydraulic conductivity [cm / h] in TopSoil

"""Physical distance to chalk/groundwater"""
chalk_depth = 150 #cm - Depth below surface that chalk occur
GW_dist = 1000 #cm - Distance to the groundwater

"""Properties in the chalk"""
theta_s_m = 0.395  # Water content at saturation, limestone matrix
theta_s_f = 0.02  # Water content at saturation, limestone fracture
b_chalk = 2  # Guess Per til Sofus' gruppe: "20-25, måske 15".
Ks_m = 6.77e-3 # cm/h - Slowest determined in Darcy experiment, (4.167 * 10^(-3) cm/h in Mathias et. al 2006)

"""Model parameters, discretization, simulation duration, step save etc"""
dt = 1  # Time step [h]
dz = 10  # Vertical cell size [cm]
t_n = 24*365*30  # Last time step [dt]

stepsave = 24  # save every nth step.

OG_values = 1
Preci = 1 #Precipitation factor
"""--------------------------------------------------"""
n_cells = GW_dist//dz  # Number of cells [-]
dtn = []
dtn2 = []
checks = []
Rt = []
vTS = []
ET_overloader = []
MB_water = []  # The transition from Topsoil to the limestone
# Water content at depth of suction cell in Gistrup (1 m) +/- one cell
PM_range = np.zeros((n_cells, 1))
# The decreasing profile of ET0, linear decrease in 50 cm depth.
PM_decrease = np.arange(1, -dz/ 50, -dz/50).reshape(int(1/(dz/50)+1), 1)
PM_decrease[-1] = 0
if ForestType == 'Pine':
    Root_depth_0 = 10  # cm
    PM_range[0:Root_depth_0//dz] = (np.ones((Root_depth_0//dz, 1)))
    PM_range[Root_depth_0//dz-1:Root_depth_0 //
        dz+len(PM_decrease)-1] = PM_decrease
    Kc = np.ones(12)*1.4
    Kc[3:9] = 1.5
    Root_growth = 5  # cm/month
    max_root = 50  # cm #Max potential of root, the following 50 cm is the decreasing Penman
elif ForestType == 'BLT':
    Root_depth_0 = 10  # cm
    PM_range[0:Root_depth_0//dz] = (np.ones((Root_depth_0//dz, 1)))
    PM_range[Root_depth_0//dz-1:Root_depth_0 //
        dz+len(PM_decrease)-1] = PM_decrease
    Kc = np.ones(12)*0.85
    Kc[3:9] = 1.05
    Root_growth = 4  # cm/month
    max_root = 100  # cm #Max potential of root, the following 50 is the decreasing Penman
else:
    Root_depth_0 = 50  # cm
    PM_range[0:Root_depth_0//dz] = (np.ones((Root_depth_0//dz, 1)))
    # PM_range[Root_depth//dz-1:Root_depth//dz+len(PM_decrease)-1] = PM_decrease
    Kc = np.ones(12)*1.05
    Root_growth = 0  # cm/month
    PM_decrease = np.zeros((6, 1))
    PM_decrease[0] = 1
    max_root = 50 #cm


# Constants

theta_wp = 0.048944  # Water content at Wilting point (Pavel 2020)

# TODO: Make soil parameters into vectors!
psi_entry = -1.87  # Air-entry pressure [cm H2O], (Pavel 2020)
b = 4.18  # 4.18 Campbells, b [-] (Pavel 2020)
# From upper Campbell interval to Wilting point (pF 3.5 to pF 4.2)
psi_interval = np.linspace(np.log10(15000), np.log10(3162), 100)
theta_interval = np.linspace(theta_wp, 0.0794, 100)
slope = (psi_interval[0]-psi_interval[-1]) / \
         (theta_interval[0]-theta_interval[-1])
skaering = psi_interval[0] - slope*theta_interval[0]



Ks_f = 4.167 # cm/h - Fastest determined in Darcy experiment (In Mathias et. al 2006 Ks_f = 41.67 cm/h)
psi_entry_chalk = -3000 # [cm H2O] Mathias et. al 2006, Per til Sofus' gruppe: "50-100 cm"

# Depth of the centroids of the cells [cm]
depth = np.array(np.arange(dz/2, dz*n_cells, dz))
topreservoir = (regn[0]-Interception[0])*dt/dz
MB = []
tr2 = []
# Initial conditions
# Water content at the first time step [cm^3 H2O / cm^3 soil] #Assumed start at pF1.5
theta_initial = 0.24
theta_fc = theta_s*(psi_entry/-100)**(1/b)  # Water content at Field capacity [cm^3 H2O / cm^3 soil ]
theta_i_m = 0.35  # Water content at first time step [cm³ H2O / cm³ ]
theta_i_f = np.zeros(n_cells)
# Soil-water content at the first time step [cm^3 H2O / cm^3 soil]
theta = np.ones(n_cells)*theta_initial
theta[chalk_depth//dz:] = theta_i_m
theta_i = np.ones(n_cells)
theta_i[0:chalk_depth//dz] = theta_initial
theta_i[chalk_depth//dz:] = theta_i_m
# Soil-water content at saturation levels [cm^3 H2O / cm^3 soil]
theta_sat = np.ones((1, n_cells))[0] * theta_s
theta_f = np.zeros(n_cells)
Nitrate = np.zeros(n_cells)
Nitrate[0] = 170  # mg/cm^3
N_save = np.copy(Nitrate)
theta_temp = np.copy(theta)
# Hydraulic conductivity in nodal points
# Hydraulic conductivity in nodal points [cm / h]
K = np.zeros((1, n_cells))[0]
K_f = np.zeros(n_cells)  # Hydraulic conductivity chalk fracture
# Saturated hydraulic conductivity in nodal points [cm / h]
K_sat = np.ones((1, n_cells))[0] * K_s
# Slope and y-intersects in nodal points
# Local y-intersect for the nodal points [cm / h]
K_L = np.zeros((1, n_cells))[0]
alpha_L = np.zeros((1, n_cells))[0]  # Local slope of the line intersecting [-]

# Slope and y-intersects in control surfaces
# Local y-intersect on the interface between cell i and i+1 [cm / h]
K_int = np.zeros((1, n_cells))[0]
# Local slope on the interface between cell i and i+1 [-]
alpha_int = np.zeros((1, n_cells))[0]

# Suction
psi = np.zeros((1, n_cells))[0]  # Soil-water potential [cm H2O]

# Velocity
v_int = np.zeros((1, n_cells))[0]  # Velocities in each nodal point [cm/h]
v_f = np.zeros(n_cells)  # velocity of chalk fracture
# v_upper = np.zeros(t_n//stepsave) # Upper boundary condition [cm/h]
# v_upper[1] = 0.1
# List to store time steps
theta_list = np.zeros((t_n//stepsave, n_cells))
theta_list[0] = theta
N_list = np.zeros((t_n//stepsave, n_cells))
K_list = np.zeros((t_n//stepsave, n_cells))
v_list = np.zeros((t_n//stepsave, n_cells))
K_L_list = np.zeros((t_n//stepsave, n_cells))
psi_list = np.zeros((t_n//stepsave, n_cells))
ET_0 = np.zeros(n_cells)
ET_Overload = np.zeros(n_cells)

ET_0_list = np.zeros((t_n//stepsave, n_cells))
K_f_list = np.zeros((t_n//stepsave, n_cells))
v_f_list = np.zeros((t_n//stepsave, n_cells))
theta_f_list = np.zeros((t_n//stepsave, n_cells))

"""Massbalance"""
P = 0  # Precipitation
ET = 0  # Evapotranspiration
I = 0  # Interception
Ud_f = 0  # Out from Fracture
Ud = 0  # Out from Matrix
N_start = sum(Nitrate*theta)*dz
N_ind = 0  # Nitrate into the system
N_out = 0


# Calculate variables
for t in range(0, t_n):
    # Calculation of variables in nodal points
    for i in range(n_cells):
        # Suction
        if i <= (chalk_depth-dz)//dz:
            if theta[i] > theta_interval[-1]:
                psi[i] = (psi_entry * (theta[i] / theta_s)**(-b))  # Campbell
            else:
                psi[i] = -10**(slope*theta_interval[i] + skaering)  # Linear
            # Hydraulic conductivity
            # Eq. (3.5) (Loll & Moldrup)
            K[i] = K_sat[i] * (psi_entry / psi[i])**(2 + b / 3)

            # Slope and y-intersect
            alpha_L[i] = (2 + 3/b) / (-psi[i])

            # Hydraulic conductivity
            K_L[i] = K[i] * np.exp(- alpha_L[i] * psi[i])

        else:
            # psi[i] = psi_entry_chalk * (theta[i] / theta_s_m) ** (-1/lambda_chalk)
            K[i] = Ks_m * (theta[i] / theta_s_m)**(2*b_chalk + 3)
            # K[i] = Ks_m * (psi_entry_chalk/psi[i])**(lambda_chalk*eta_chalk)

    # Calculation of variables at the control surfaces
    for i in range(0, n_cells):
        # Calculate values for interfaces between nodal points (control surfaces)
        if i == n_cells-1 or i == (chalk_depth-dz)//dz:
            # Lower boundary condition
            v_int[i] = K[i]
        elif i >= chalk_depth//dz:
            v_int[i] = K[i]
        else:
            # Slope
            alpha_int[i] = (alpha_L[i] + alpha_L[i+1]) / 2

            # Hydraulic conductivity
            K_int[i] = (K_L[i] + K_L[i+1]) / 2

            # Velocity
            K_n1 = K_int[i] * np.exp(alpha_int[i] * psi[i])
            K_n2 = K_int[i] * np.exp(alpha_int[i] * psi[i+1])

            v_int[i] = - (K_n2 - K_n1) / (np.exp(alpha_int[i] * dz) -
                          1) + K_int[i] * np.exp(alpha_int[i] * psi[i])

    # We use the water content as a initial conditions, so only calculate the
    # water content for the rest of the time steps

    if t != 0:
        # Calculate the new water content in nodal points
        for i in range(0, n_cells):
            
           # """If theta[i] == theta_wp:
            #   theta[i] = theta[i] + Σ - 0, else theta[i] = theta[i] + Σ - ET"""
            # In the first box (i =) 0, we need to specify the velocity from the
            # upper boundary condition (v_upper)
            #if theta[i] == theta_wp:  # if theta is at wilting point
                # Temporary save of ET at given time
                #ET_temp = ET0.iloc[t % len(ET0)]
                # New value of ET, since the soil is at wilting point.
                #ET0.iloc[t % len(ET0)] = 0
                #v_int[i] = 0
               # Checker = True
            #else:
               # Checker = False
            PM_range[0:Root_depth//dz] = (np.ones((Root_depth//dz, 1)))
            PM_range[Root_depth//dz-1:Root_depth //
                dz+len(PM_decrease)-1] = PM_decrease
            ET_0[i] = Kc[mon[t % len(regn)]] * ET0.iloc[t % len(ET0)]*PM_range[i]/(PM_range.sum())*min((theta[i]-theta_wp)/(theta_fc-theta_wp), 1)
            if i == 0:
                topreservoir2 = topreservoir
                topreservoir = topreservoir - (min(topreservoir*dz/dt*0.5, 0.5) -
                     regn[t % len(regn)]+Interception[t % len(Interception)])*dt/dz

                # topreservoir = max((topreservoir - (K_s*topreservoir + Kc[mon[t%len(regn)]] * ET0.iloc[t%len(ET0)]*PM_range[i,0]/(PM_range.sum()+1)*min((theta[i]-theta_wp)/(theta_initial-theta_wp), 1) - regn[t%len(regn)]) * dt / dz), 0)
                tr2.append(topreservoir2)
                New_temp_theta = np.copy(theta[i] - (v_int[i] + ET_0[i] -
                               min(topreservoir2*dz/dt*0.5, 0.5)) * dt / dz)
                
                theta[i] = max((theta[i] - (v_int[i] + ET_0[i] -
                               min(topreservoir2*dz/dt*0.5, 0.5)) * dt / dz), theta_wp)
                if v_int[i] < 0:
                    Nitrate[i] = max((-(v_int[i]*N_save[i+1])*dt/dz+theta_temp[i] * N_save[i])/theta[i], 1e-20)
                else:
                    Nitrate[i] = max((-(v_int[i]*N_save[i])*dt/dz+theta_temp[i] * N_save[i])/theta[i], 1e-20)
               # ET_Overload[i] = min(New_temp_theta-theta_wp, 0)


            # For the rest of the boxes, we use the calculated velocities at the
            # control surfaces
                # [t%len(regn)])
           # *topreservoir2  / dz
            # *topreservoir/dz
            elif i >= chalk_depth//dz:
                New_temp_theta = theta[i] - (v_int[i] - v_int[i-1])* dt/dz + theta_f[i-1]
                theta[i]=min(theta[i] - (v_int[i] - v_int[i-1])
                             * dt/dz + theta_f[i-1], theta_s_m)
                
                # Since theta_f[t-1,i] is fully drained, theta_f[i] - theta_f[i] is neglected
                theta_f[i]=max(New_temp_theta-theta_s_m, 0)
                v_f[i] = theta_f[i]*dz/dt 
                Nitrate[i]=(-(v_int[i]*Nitrate[i]-v_int[i-1]*N_save[i-1])*dt/dz+theta_temp[i]*Nitrate[i])/theta[i]
                
            else:
                New_temp_theta = (theta[i] - (v_int[i] - v_int[i-1]+ ET_0[i]) * dt / dz)
                theta[i]=max(
                    (theta[i] - (v_int[i] - v_int[i-1] + ET_0[i]) * dt / dz), theta_wp)
                #ET_Overload[i] = min(New_temp_theta-theta_wp, 0)
                if v_int[i] < 0 and v_int[i-1] < 0:
                    Nitrate[i]= max((-(v_int[i]*N_save[i+1]-v_int[i-1]*N_save[i])*dt/dz+theta_temp[i]*N_save[i])/theta[i], 1e-20)
                elif v_int[i] < 0:
                    Nitrate[i]= max((-(v_int[i]*N_save[i+1]-v_int[i-1]*N_save[i-1])*dt/dz+theta_temp[i]*N_save[i])/theta[i], 1e-20)
                elif v_int[i-1] < 0:
                    Nitrate[i]= max((-(v_int[i]*N_save[i]-v_int[i-1]*N_save[i])*dt/dz+theta_temp[i]*N_save[i])/theta[i], 1e-20)
                else:
                    Nitrate[i]= max((-(v_int[i]*N_save[i]-v_int[i-1]*N_save[i-1])*dt/dz+theta_temp[i]*N_save[i])/theta[i], 1e-20)
            #if Checker == True:  # If ET0 was changed to 0, change it back to the temporary saved value.
             #   ET0.iloc[t % len(ET0)]=ET_temp
                
        # Temporary save of theta, so secure the new concentration keeps balance in the mass.
        ET_overloader.append(ET_Overload)
        theta_temp=np.copy(theta)
        N_save=np.copy(Nitrate)
        MB.append(sum(Nitrate*theta)*dz)
        
        if t % stepsave == 0:
            theta_list[t//stepsave]=theta

            ET_0_list[t//stepsave]=ET_0
            theta_f_list[t//stepsave]=theta_f
        dtn.append(min(theta))
        dtn2.append(max(theta_f))
        # if i < n_cells:
           # dtn[i] = ((theta_list[i+1]-theta_list[i])/(K_list[i+1]-K_list[i]))*(np.exp(alpha_int[i]*dz-1))/(np.exp(alpha_int[i]*dz+1))*dz

    if t < len(mon)-100:
        Root_depth=min(Root_depth_0 + Root_growth * \
                       mon[t % len(mon)] + yr[t % len(yr)] * Root_growth * 12, max_root)
        Rt.append(Root_depth)


    # Store variables in nested arrays
    if t % stepsave == 0:
        K_list[t//stepsave]=K
        v_list[t//stepsave]=v_int
        K_L_list[t//stepsave]=K_L
        psi_list[t//stepsave]=psi
        N_list[t//stepsave]=Nitrate
        v_f_list[t//stepsave]=v_f
        theta_f_list[t//stepsave] = theta_f 


    vTS.append(v_int[(chalk_depth-dz)//dz])
    P=P + regn[t % len(regn)] * dt
    ET=ET + (sum(ET_0)) * dt
    Ud=Ud + v_int[-1] * dt
    I=I + Interception[t % len(Interception)] * dt
    Ud_f=Ud_f + theta_f[-1]*dz
    N_out=N_out + v_int[-1] * dt * Nitrate[-1]
    PM.append(sum(ET_0))
    MB_check.append(P-(ET+Ud+I+Ud_f+(sum(theta)-sum(theta_i))*dz-topreservoir*dz))

    if t == t_n-1:
        M=(sum(theta)-sum(theta_i))*dz  # Moisture topsoil + limestone matrix
        F=(sum(theta_f)-sum(theta_i_f))*dz  # Moisture fracture
        N_akk=sum(Nitrate*theta)*dz-N_start

print(datetime.now()-begin_time)
print(min(dtn))
print(max(dtn2))
print(' Waterbalance: \n Precipitation: ' + str(P) + '\n Interception: ' + str(I) + '\n Evapotranspiration: ' + str(ET) + '\n Out matrix: '
      + str(Ud) + '\n Out Fracture: ' + str(Ud_f) + \
            '\n Δ(Soil + Matrix moisture): '
      + str(M) + '\n ΔFracture moisture: ' + str(F))
print('P-(I+ET+ΔS+Out) = ' + str(P-(I+ET+M+F+Ud+Ud_f)))

print('Nitrate massbalance:' + ' \n In: ' + str(N_ind) + '\n Out: ' + str(N_out) + \
      '\n Accumulated: ' + str(N_akk) + '\n Sum: ' + str(N_ind - (N_out + N_akk)))
MB = np.array(MB)
vTS = np.array(vTS)
np.savez(ForestType+'-'+ParameterChange+'-'+str(globals()[ParameterChange])+'.npz', N_liste = N_list, ET0_liste = ET_0_list, K_liste = K_list, theta_liste = theta_list, theta_f_liste = theta_f_list, v_frac_liste = v_f_list, v_liste = v_list, MB_nitrat = MB, v_TS = vTS, Precipitation = P, Interception = I, Out_matrix = Ud, Storage = M, ET = ET, Out_frac = Ud_f, Storage_Frac = F)
print((pd.DataFrame(N_list).min() < 0).sum())
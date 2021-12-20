# Afforestation Project
The afforestation project is a semester project on the first semester of the master in 'Water and Environmental Engineering' at the institute of Environemntal Built at Aalborg University.

The numerical model is made to determine the effect of afforestation in Gistrup, Denmark and how the afforestation effects the groundwater recharge, combined with the nitrate leaching, and is based on Python 3.8.8.

The file 'Numerical_Afforestation_Model' is the main file, that requires some input from other files, stored as csv-files.
The csv-files are based on:
  1) Precipitation, measured at Oesterport, Aalborg in Denmark, sliced to a 4-year period (2017-2020)
  2) Measured weather data from the Danish Meteorological Institute for the 4-year period (2017-2020). The measured data are used to calculate the Penman-Monteith equation.
  3) Interception, based on a forecast from the Leaf Area Index (LAI), and then adjusted to the precipitation data. 


The precipitation and calculation of the Penman-Monteith are stored in 'DMISamletDataFrame.csv', while the output of the Penman-Monteith equation is stored in the 'Penman_"vegetation".csv', where "vegetation" corresponds to the used vegetation type. Grass for Grassland, BLT for a Deciduous forest and Pine for a Coniferous forest.
The csv files for Interception is 'Interception_"vegetation".csv'.

The '5052Oesterport_alle_regn_1979-2021.km2' file is the measured rain at Oesterport for the period 1979-2021.

The model calculates and outputs NumPy arrays that are used to visualize the results.

In order to produce the csv-files that are required in the numerical model, the 'Interception.py' and 'Penman.py' are run.


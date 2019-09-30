#!/usr/bin/env python
# coding: utf-8


import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii

import rasterio
from rasterio import plot

get_ipython().run_line_magic('matplotlib', 'inline')
import earthpy as et
import earthpy.plot as ep
import copy
import pandas as pd
import xarray as xr

#Path to folder with our data
path_to_temp_mean_folder = "./Huglin_Index_Data/CMIP5_Downscaled_GDF_ds01_1-16/cmip5dc_global_10min_rcp2_6_2030s_asc/bcc_csm1_1_m_rcp2_6_2030s_tmean_10min_r1i1p1_no_tile_asc/"
path_to_temp_max_folder = "./Huglin_Index_Data/CMIP5_Downscaled_GDF_ds01_1-16/cmip5dc_global_10min_rcp2_6_2030s_asc/bcc_csm1_1_m_rcp2_6_2030s_tmax_10min_r1i1p1_no_tile_asc/"

no_data_value = -9999


#Creating two lists of mean and max temperatures of the 6 months we examine.
temp_mean_list = []
temp_max_list = []
#Months index i, April to September (4 to 9)
for i in range(4,10):
    #Getting the data without the header lines and as float type
    #they are stored as int and multiplied by 10 to preserve the first decimal 
    t_mean_array = (np.loadtxt(path_to_temp_mean_folder + "tmean_" + str(i) + ".asc", skiprows=6)).astype('float64')
    #Data values divided by 10 to get their real value, no data values are ignored
    t_mean_array[t_mean_array != no_data_value] = t_mean_array[t_mean_array != no_data_value]/10

    temp_mean_list.append(t_mean_array)
    
    
    
    t_max_array = (np.loadtxt(path_to_temp_max_folder + "tmax_" + str(i) + ".asc", skiprows=6)).astype('float64')
    t_max_array[t_max_array != no_data_value] = t_max_array[t_max_array != no_data_value]/10
    
    temp_max_list.append(t_max_array)



#Calculating party Huglin Index  
huglin_index = np.zeros(temp_mean_list[0].shape, dtype=float)
for i in range(len(temp_mean_list)):
    #Checking which months have 30 or 31 days and apply the equation
    if i==4 or i ==6 or i==9:
        huglin_index = huglin_index + 30*(temp_mean_list[i] + temp_max_list[i] - 20)/2
    else:
        huglin_index = huglin_index  + 31*(temp_mean_list[i] + temp_max_list[i] - 20)/2




#Checking properties of the result array
print(huglin_index.shape)



huglin_min_value = np.amin(huglin_index)
print(huglin_min_value)

#Correcting the values of the no data values to the global value -9999
huglin_index[huglin_index < no_data_value] = no_data_value

#Getting mim/max values ignoring no data values
huglin_max_value = np.amax(huglin_index[huglin_index != no_data_value])
huglin_min_value = np.amin(huglin_index[huglin_index != no_data_value])
print(huglin_min_value)
print(huglin_max_value)


#Visualizing the result so far
ep.plot_bands(huglin_index,cmap='PiYG',scale=False, vmin=huglin_min_value,vmax=huglin_max_value)



#Creating the latitude, longtitude 
x,y = huglin_index.shape
lat = np.linspace(-90, 90, x)
lon = np.linspace(-180, 180, y)

#Creating dataset using the latitude/longtitude
huglin_dataset = xr.DataArray(data=huglin_index, dims=["lat", "lon"], coords=[lat,lon])
print(huglin_dataset)

#Getting an entry using lat lon
test = huglin_dataset.sel(lat=-4.2,lon=-9.7,method='nearest')   
print(test.data)




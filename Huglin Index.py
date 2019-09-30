#!/usr/bin/env python
# coding: utf-8

# In[1]:


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

path_to_temp_mean_folder = "../Huglin_Index_Data/CMIP5_Downscaled_GDF_ds01_1-16/cmip5dc_global_10min_rcp2_6_2030s_asc/bcc_csm1_1_m_rcp2_6_2030s_tmean_10min_r1i1p1_no_tile_asc/"
path_to_temp_max_folder = "../Huglin_Index_Data/CMIP5_Downscaled_GDF_ds01_1-16/cmip5dc_global_10min_rcp2_6_2030s_asc/bcc_csm1_1_m_rcp2_6_2030s_tmax_10min_r1i1p1_no_tile_asc/"

no_data_value = -9999


# In[2]:


temp_mean_list = []
temp_max_list = []
for i in range(4,10):
    t_mean_array = (np.loadtxt(path_to_temp_mean_folder + "tmean_" + str(i) + ".asc", skiprows=6)).astype('float64')
    t_mean_array[t_mean_array != no_data_value] = t_mean_array[t_mean_array != no_data_value]/10

    temp_mean_list.append(t_mean_array)
    
    
    
    t_max_array = (np.loadtxt(path_to_temp_max_folder + "tmax_" + str(i) + ".asc", skiprows=6)).astype('float64')
    t_max_array[t_max_array != no_data_value] = t_max_array[t_max_array != no_data_value]/10
    
    temp_max_list.append(t_max_array)



# In[4]:


huglin_index_list = []
huglin_index = np.zeros(temp_mean_list[0].shape, dtype=float)
for i in range(len(temp_mean_list)):
    if i==4 or i ==6 or i==9:
        huglin_index_list.append(30*(temp_mean_list[i] + temp_max_list[i] - 20)/2)
    else:
        huglin_index_list.append(31*(temp_mean_list[i] + temp_max_list[i] - 20)/2)
        
    huglin_index = huglin_index + huglin_index_list[i]


# In[41]:


print(huglin_index.shape)
huglin_max_value = np.amax(huglin_index[huglin_index != np.amax(huglin_index)])


huglin_min_value = np.amin(huglin_index)
print(huglin_min_value)
huglin_index[huglin_index < no_data_value] = no_data_value
print(huglin_index)


huglin_min_value = np.amin(huglin_index[huglin_index != np.amin(huglin_index)])
print(huglin_min_value)
print(huglin_max_value)


# In[7]:


huglin_index_min_value = np.amin(huglin_index[huglin_index != no_data_value])
huglin_index_max_value = np.amax(huglin_index[huglin_index != no_data_value])
print(huglin_index_min_value)
print(huglin_index_max_value)


# In[8]:


ep.plot_bands(huglin_index,cmap='PiYG',scale=False, vmin=huglin_min_value,vmax=huglin_max_value)


# In[17]:


x,y = huglin_index.shape
lat = np.linspace(-90, 90, x)
lon = np.linspace(-180, 180, y)

huglin_dataset = xr.DataArray(data=huglin_index, dims=["lat", "lon"], coords=[lat,lon])
print(huglin_dataset)


# In[35]:


test = huglin_dataset.sel(lat=-4.2,lon=-9.7,method='nearest')    #[lat: -90, lon: -180]


# In[37]:


print(test.data)


# In[38]:


test.data = 0
print(test.data)


#!/usr/bin/env python
# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii

import rasterio
from rasterio import plot

from osgeo import ogr, osr, gdal

#get_ipython().run_line_magic('matplotlib', 'inline')
import earthpy as et
import earthpy.plot as ep
import copy
import pandas as pd
import xarray as xr

path_to_temp_mean_folder = "./Huglin_Index_Data/CMIP5_Downscaled_GDF_ds01_1-16/cmip5dc_global_10min_rcp2_6_2030s_asc/bcc_csm1_1_m_rcp2_6_2030s_tmean_10min_r1i1p1_no_tile_asc/"
path_to_temp_max_folder = "./Huglin_Index_Data/CMIP5_Downscaled_GDF_ds01_1-16/cmip5dc_global_10min_rcp2_6_2030s_asc/bcc_csm1_1_m_rcp2_6_2030s_tmax_10min_r1i1p1_no_tile_asc/"

no_data_value = -9999



class Huglin_Index():
    def __init__(self, hemisphere,no_data_value):
        self.hemisphere = hemisphere
        self.no_data_value = no_data_value
        self.temp_mean_list = []
        self.temp_max_list = []
        
        
    def get_data(self,path_to_temp_mean_folder,path_to_temp_max_folder):
        
        if self.hemisphere=="N":
            
            for i in range(4,10):
                #Getting the data without the header lines and as float type
                #they are stored as int and multiplied by 10 to preserve the first decimal 
                t_mean_array = (np.loadtxt(path_to_temp_mean_folder + "tmean_" + str(i) + ".asc", skiprows=6)).astype('float64')
                #Data values divided by 10 to get their real value, no data values are ignored
                t_mean_array[t_mean_array != self.no_data_value] = t_mean_array[t_mean_array != self.no_data_value]/10

                self.temp_mean_list.append(t_mean_array)

                t_max_array = (np.loadtxt(path_to_temp_max_folder + "tmax_" + str(i) + ".asc", skiprows=6)).astype('float64')
                t_max_array[t_max_array != self.no_data_value] = t_max_array[t_max_array != self.no_data_value]/10

                self.temp_max_list.append(t_max_array)
        else:
            
            for i in range(1,4):
                t_mean_array = (np.loadtxt(path_to_temp_mean_folder + "tmean_" + str(i) + ".asc", skiprows=6)).astype('float64')
                t_mean_array[t_mean_array != self.no_data_value] = t_mean_array[t_mean_array != self.no_data_value]/10

                self.temp_mean_list.append(t_mean_array)

                t_max_array = (np.loadtxt(path_to_temp_max_folder + "tmax_" + str(i) + ".asc", skiprows=6)).astype('float64')
                t_max_array[t_max_array != self.no_data_value] = t_max_array[t_max_array != self.no_data_value]/10

                self.temp_max_list.append(t_max_array) 
                
            for i in range(10,13):
                t_mean_array = (np.loadtxt(path_to_temp_mean_folder + "tmean_" + str(i) + ".asc", skiprows=6)).astype('float64')
                t_mean_array[t_mean_array != self.no_data_value] = t_mean_array[t_mean_array != self.no_data_value]/10

                self.temp_mean_list.append(t_mean_array)

                t_max_array = (np.loadtxt(path_to_temp_max_folder + "tmax_" + str(i) + ".asc", skiprows=6)).astype('float64')
                t_max_array[t_max_array != self.no_data_value] = t_max_array[t_max_array != self.no_data_value]/10

                self.temp_max_list.append(t_max_array)
        return
                
    def calculate_index(self):
        self.huglin_index = np.zeros(self.temp_mean_list[0].shape, dtype=float)
        if self.hemisphere=="N":
            month = 4
            for i in range(len(self.temp_mean_list)):
                if month==4 or month==6 or month==9:
                    huglin_index_var = 30*(self.temp_mean_list[i] + self.temp_max_list[i] - 20)/2
                else:
                    huglin_index_var = 31*(self.temp_mean_list[i] + self.temp_max_list[i] - 20)/2

                self.huglin_index = self.huglin_index  + huglin_index_var
                month += 1
        else:
            month = 1
            for i in range(len(self.temp_mean_list)//2):
                if month==2:
                    huglin_index_var = 28*(self.temp_mean_list[i] + self.temp_max_list[i] - 20)/2
                else:
                    huglin_index_var = 31*(self.temp_mean_list[i] + self.temp_max_list[i] - 20)/2

                self.huglin_index = self.huglin_index  + huglin_index_var
                month += 1
                
            month = 10
            for i in range(len(self.temp_mean_list)//2):
                if month==11:
                    huglin_index_var = 30*(self.temp_mean_list[i] + self.temp_max_list[i] - 20)/2
                else:
                    huglin_index_var = 31*(self.temp_mean_list[i] + self.temp_max_list[i] - 20)/2

                self.huglin_index = self.huglin_index  + huglin_index_var
                month += 1
                
        
                
        x,y = self.huglin_index.shape
        self.lat = np.linspace(-180, 180, x)
        self.lon = np.linspace(-90, 90, y)

        self.huglin_dataset = xr.DataArray(data=self.huglin_index, dims=["lat", "lon"], coords=[self.lat,self.lon])
        
        lat_array_indexes= []
        if self.hemisphere == "N":
            lat0 = 40
            for i in range(0,11,2):
                tmp = float(self.huglin_dataset.sel(lat=lat0 + i,lon=-90,method='nearest').lat.data )
                lat_array_indexes.append(int(np.where(self.lat==tmp)[0]))
        else:
            lat0 = -40
            for i in range(0,-11,-2):
                tmp = float(self.huglin_dataset.sel(lat=lat0 + i,lon=-90,method='nearest').lat.data )
                lat_array_indexes.append(int(np.where(self.lat==tmp)[0]))

        
        
        k=[1.02,1.03,1.04,1.05,1.06]    
        for i in range(len(lat_array_indexes)-1):
            for j in range(lat_array_indexes[i+1] - lat_array_indexes[i]):
                self.huglin_dataset[lat_array_indexes[i] + j] = self.huglin_dataset[lat_array_indexes[i] + j] * k[i]
                #pointer_to_row_latitude = self.huglin_dataset[lat_array_indexes[i] + j]
                #pointer_to_row_latitude[self.huglin_dataset[lat_array_indexes[i] + j] < no_data_value] = no_data_value
        
        
        self.huglin_max = np.amax(self.huglin_dataset.data)
        self.huglin_dataset.data[self.huglin_dataset < no_data_value] = no_data_value

        self.huglin_min = np.amin(self.huglin_dataset.data[self.huglin_dataset != np.amin(self.huglin_dataset)])

        
        return self.huglin_dataset
                
                
                
    def show_map(self):
        if self.hemisphere == "N":
            print("Huglin Index of North Hemisphere")
        else:
            print("Huglin Index of North Hemisphere")
        print("Max value of Huglin Index: " + str(self.huglin_max))
        print("Min value of Huglin Index: " + str(self.huglin_min))
        ep.plot_bands(self.huglin_dataset,cmap='PiYG',scale=False, vmin=self.huglin_min,vmax=self.huglin_max)
        
        return
    
    def return_huglin_index_from_lat_lon(self, lat, lon):
        return self.huglin_dataset.sel(lat=lat,lon=lon,method='nearest')
    
    
    def CreateGeoTiff(self, outRaster):
        data = []
        data.append(self.huglin_index)
        data.extend(self.temp_mean_list)
        data.extend(self.temp_max_list)
        data = np.array(data)
        
        driver = gdal.GetDriverByName('GTiff')
        no_bands, width, height = data.shape
        DataSet = driver.Create(outRaster, height, width, no_bands, gdal.GDT_Float64)
        DataSet.SetGeoTransform([-180,0.1666666666667,0,-60,0,0.1666666666667])

        for i, image in enumerate(data, 1):
            DataSet.GetRasterBand(i).WriteArray(image)
        DataSet = None
       
        return



north_hem_huglin = Huglin_Index("N",no_data_value)
south_hem_huglin = Huglin_Index("S",no_data_value)

  


north_hem_huglin.get_data(path_to_temp_mean_folder,path_to_temp_max_folder)
south_hem_huglin.get_data(path_to_temp_mean_folder,path_to_temp_max_folder)



north_hem_huglin_index_object = north_hem_huglin.calculate_index()
south_hem_huglin_index_object = south_hem_huglin.calculate_index()

north_hem_huglin_index = north_hem_huglin_index_object.data
south_hem_huglin_index = south_hem_huglin_index_object.data



north_hem_huglin.show_map()
north_hem_huglin.return_huglin_index_from_lat_lon(13,-40)

south_hem_huglin.show_map()
south_hem_huglin.return_huglin_index_from_lat_lon(13,-40)

north_hem_huglin.CreateGeoTiff("./Huglin_Index_Mean_Max_Temps.tiff")

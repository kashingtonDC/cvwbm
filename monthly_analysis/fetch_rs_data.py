import os
import ee
import datetime
import time

import geopandas as gp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from dateutil.relativedelta import relativedelta

ee.Initialize()

'''
Set the init params:
1) Area: Read from Shapefile
2) Years to analyze: Supplied
3) Scale (optional): Resolution at which to perform analysis
4) Satellites / datasets: Queried from Earth Engine in the datastructure below: 
'''

data = {}

###################
##### ET data #####
###################

data['modis'] = [ee.ImageCollection('MODIS/006/MOD16A2'), "ET", 0.1]
data['gldas'] = [ee.ImageCollection('NASA/GLDAS/V021/NOAH/G025/T3H'), 'Evap_tavg', 86400]
data['nldas'] = [ee.ImageCollection('NASA/NLDAS/FORA0125_H002'), 'potential_evaporation', 1]
data['terra'] = [ee.ImageCollection('IDAHO_EPSCOR/TERRACLIMATE'), "aet", 0.1]
data['gmet_etr'] = [ee.ImageCollection('IDAHO_EPSCOR/GRIDMET'), "etr", 1]
data['gmet_eto'] = [ee.ImageCollection('IDAHO_EPSCOR/GRIDMET'), "eto", 1]
data['fldas'] = [ee.ImageCollection('NASA/FLDAS/NOAH01/C/GL/M/V001'), "Evap_tavg", 86400*24 ]
data['cas'] = [ee.ImageCollection('CAS/IGSNRR/PML/V2'), "ET_water", 1]

###################
##### P data ######
###################

data['trmm']  =  [ee.ImageCollection('TRMM/3B43V7'), "precipitation", 720]
data['prism'] = [ee.ImageCollection("OREGONSTATE/PRISM/AN81m"), "ppt", 1]
data['chirps'] = [ee.ImageCollection('UCSB-CHG/CHIRPS/PENTAD'), "precipitation", 1]
data['persia'] = [ee.ImageCollection("NOAA/PERSIANN-CDR"), "precipitation", 1]
data['dmet'] = [ee.ImageCollection('NASA/ORNL/DAYMET_V3'), "prcp", 1]

###################
##### SWE data #####
####################
data['swe'] = [ee.ImageCollection('NASA/FLDAS/NOAH01/C/GL/M/V001'), "SWE_inst", 1 ]

####################
##### DS data ######
####################
data['grace'] = [ee.ImageCollection('NASA/GRACE/MASS_GRIDS/LAND'), "lwe_thickness_jpl", 1]

####################
##### R data #######
####################
data['r'] = [ee.ImageCollection('IDAHO_EPSCOR/TERRACLIMATE'), "ro", 1]
data['runoff'] = [ee.ImageCollection("NASA/FLDAS/NOAH01/C/GL/M/V001"), "Qs_tavg", 86400*24]

#####################
##### SM data #######
#####################
data['sm'] = [ee.ImageCollection('IDAHO_EPSCOR/TERRACLIMATE'), "soil", 0.1]
data['sm1'] = [ee.ImageCollection("NASA/FLDAS/NOAH01/C/GL/M/V001"), "SoilMoi00_10cm_tavg", 1]
data['sm2'] = [ee.ImageCollection("NASA/FLDAS/NOAH01/C/GL/M/V001"), "SoilMoi10_40cm_tavg", 1]
data['sm3'] = [ee.ImageCollection("NASA/FLDAS/NOAH01/C/GL/M/V001"), "SoilMoi40_100cm_tavg", 1]
data['sm4'] = [ee.ImageCollection("NASA/FLDAS/NOAH01/C/GL/M/V001"), "SoilMoi100_200cm_tavg", 1]


# Helper functions to get the relevant EE dataset over AOI. In this case, get the sum, convert to to km^3, return 

def get_data(dataset, year, month):

    col = dataset[0]
    var = dataset[1]
    scaling_factor = dataset[2]

    t = col.filter(ee.Filter.calendarRange(year, year, 'year')).filter(ee.Filter.calendarRange(month, month, 'month')).select(var).filterBounds(area).sum()
    t2 = t.multiply(1e-3).multiply(ee.Image.pixelArea()).multiply(scaling_factor).multiply(1e-9)
    # convert mm to m, multiply by pixel area (m^2), multiply by scaling factor, m^3 to km^3
    scale = t2.projection().nominalScale()
    sumdict  = t2.reduceRegion(
        reducer = ee.Reducer.sum(),
        geometry = area,
        scale = scale)
    
    result = sumdict.getInfo()[var]
    
    return result

def wrapper(dataset):
	# TODO: Error handling try/except statement here 
    monthly = []

    for year in years:
        for month in months:
            r = get_data(dataset, year, month)
            monthly.append(r)
            time.sleep(5) # in case this hits ee memory limits
    
    print("wrapper complete")
    return monthly

def read_shapefile(shapefile):

	# get the area from a shapefile
	cv = gp.read_file(shapefile)
	cv2 = cv.geometry.simplify(0.1)
	lls = cv2.geometry.iloc[0]
	x,y = lls.exterior.coords.xy
	coords = [list(zip(x,y))]
	area = ee.Geometry.Polygon(coords)
	sa_km2 = round(cv2.area[0] * 10000, 2)

	return area


def main(shapefile):

	global area 
	area = read_shapefile(shapefile)

	# TODO: Make the years and datasets be read from a parameter file? 

	years = range(2001, 2019)
	months = range(1,13)
	start = datetime.datetime(years[0], 1, 1)
	end = datetime.datetime(years[-1]+1, 1, 1)
	dt_idx = pd.date_range(start,end, freq='M')

	# fetch the data

	# soil moisture
	sm = wrapper(data['sm'])
	# swe
	swe = wrapper(data['swe'])
	# runoff
	r = wrapper(data['r'])
	# et
	et_m = wrapper(data['modis'])
	eto_g = wrapper(data['gmet_eto'])
	etr_g = wrapper(data['gmet_etr'])
	et_f = wrapper(data['fldas'])
	et_t = wrapper(data['terra'])
	# TODO: gldas, nldas

	# precipitation
	p_p = wrapper(data['prism'])
	p_c = wrapper(data['chirps'])
	p_n = wrapper(data['persia'])
	p_t = wrapper(data['trmm'])
	p_d = wrapper(data['dmet'])

	# Make calculate temporal means when you have more than one dataset (e.g. P, ET) 
	p_arrs = [np.array(x) for x in [p_p, p_c, p_n, p_t, p_d]]
	p_stacked = np.column_stack(p_arrs)
	p = np.mean(p_stacked, axis = 1)

	et_arrs = [np.array(x) for x in [et_m, eto_g, etr_g, et_f, et_t]]
	et_stacked = np.column_stack(et_arrs)
	et = np.mean(et_stacked, axis = 1)

	# Make a df from the RS data:
	years = range(2001, 2018)
	months = range(1,13)
	start = datetime.datetime(years[0], 1, 1)
	end = datetime.datetime(years[-1]+1, 1, 1)
	dt_idx = pd.date_range(start,end, freq='M')

	rsdf = pd.DataFrame([r, sm, p_p, p_c, p_n, p_t,p_d, et_m, eto_g, etr_g, et_f, et_t, list(p),list(et)]).T

	rsdf.index = dt_idxs

	df.columns = ["runoff", "soil_moisture", "prism", "chirps", 
              "persiann" ,"trmm", "daymet", "modis", "gmet_eto", 
              "gmet_etr", "fldas", "terraclimate", 'p_mean', 'et_mean']

	df.to_csv("rs_data_monthly.csv")

if __name__ == "__main__":
	years = [x for x in range(2001, 2018)]
	months = range(1,13)
	shapefile  = "HU4_merged.shp"
	main(shapefile)
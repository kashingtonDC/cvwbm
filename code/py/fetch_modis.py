import sys
import os
import ee
import eefuncs as eef
from eefuncs import array_from_df
import numpy as np

ee.Initialize()

if __name__ == "__main__": 
	if len(sys.argv) != 3:
		# incorrect argumemnts, print usage message
		print("Usage:")
		print("  $ python3 fetch_modis.py <county_number> <start_doy>")
		sys.exit(0)

	# Print the supplied args
	print("Processing county number: "+ sys.argv[1])
	county_number = int(sys.argv[1])
	print("For day of the year: "+ sys.argv[2])
	start_doy = int(sys.argv[2])

# Functions
def filter_modis_sr(imcol, year, start_doy = start_doy):
    filtered = imcol.filter(ee.Filter.calendarRange(start_doy-5,start_doy+5, 'day_of_year')).filter(ee.Filter.calendarRange(year,year,'year')).mosaic()
    return filtered
    
def filter_modis_lc(imcol, year):
    annual_lc = imcol.filter(ee.Filter.calendarRange(year,year,'year')).select(["LC_Type1"]).median()
    return annual_lc

# Main subroutine 
def get_sr_data(col, year, start_doy, bounds):

	sr = filter_modis_sr(col, year)

	clipped_sr = ee.ImageCollection(sr.clip(bounds))

	sr_out = clipped_sr.getRegion(bounds,500,"epsg:4326").getInfo()
	t = eef.df_from_imcol(sr_out)

	# MODIS bandnames 
	bandnames = ["Nadir_Reflectance_Band1","Nadir_Reflectance_Band2","Nadir_Reflectance_Band3", "Nadir_Reflectance_Band4", "Nadir_Reflectance_Band5","Nadir_Reflectance_Band6", "Nadir_Reflectance_Band7"]

	ims = []
	for b in bandnames:
	    ims.append(array_from_df(t,b)*0.001)

	out = [np.nan_to_num(x) for x in ims]

	return out

def get_lc_data(col, year, bounds):
	lc = filter_modis_lc(col, year)
	clipped_lc = ee.ImageCollection(lc.clip(bounds))
	lc_out = clipped_lc.getRegion(bounds,500,"epsg:4326").getInfo()
	
	tlc = eef.df_from_imcol(lc_out)

	# Select the lc dataset
	landcover = array_from_df(tlc, "LC_Type1")
	# Turn ag areas into 1, others into 0
	lc = np.where(landcover == 12 , landcover, 0) 
	lc[lc>0] = 1

	out = np.nan_to_num(lc)

	return out


# main routine
def main():
	# Set the Study area (upoad kml to google fusion table, use the DocID in the ft string below)
	padded_cn = str(county_number).zfill(3)
	area = (ee.FeatureCollection('ft:1QPasan0i6O9uUlcYkjqj91D7mbnhTZCmzS4t7t_g').filter(ee.Filter().eq('id', padded_cn)))
	bounds = area.geometry().bounds()

	# Set the study years
	years = [x for x in range(2001, 2017)]
	band_names = ["b1", "b2", "b3", "b4", "b5", "b6", "b7", "mask"]

	modis_lc = ee.ImageCollection('MODIS/006/MCD12Q1')
	modis_sr = ee.ImageCollection('MCD43A4')

	for year in years:
		print("Processing: " + str(year))

		# Setup outfiles like: /modis/county_year_doy_band.txt
		moddir = os.path.join(os.getcwd(),"modis")

		if not os.path.exists(moddir):
			os.mkdir(moddir)

		filepath = os.path.join(moddir,"_".join([str(county_number),str(year), str(start_doy)]))

		outpaths = [filepath+x+".txt" for x in band_names]

		# call earth engine to get the data 
		outdata = get_sr_data(modis_sr,year, start_doy, bounds)
		outdata.append(get_lc_data(modis_lc,year, bounds))
		
		# save
		for i in range(len(outdata)):
			np.savetxt(outpaths[i], outdata[i])

		print(" === Done ===")

main()


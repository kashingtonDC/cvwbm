import sys
import os
import ee
import eefuncs as eef
from eefuncs import array_from_df, array_from_col, array_from_coords
import numpy as np

ee.Initialize()

if __name__ == "__main__": 
	if len(sys.argv) != 3:
		# incorrect argumemnts, print usage message
		print("Usage:")
		print("  $ python3 fetch_landsat.py <county_number> <start_doy>")
		sys.exit(0)

	# Print the supplied args
	print("LANDSAT")
	print("Processing county number: "+ sys.argv[1])
	county_number = int(sys.argv[1])
	print("For day of the year: "+ sys.argv[2])
	start_doy = int(sys.argv[2])

# Functions
def get_landsat(year):
	
	'''
	select the appropriate landsat based on operation years 
	'''
	
	landsats = {"L4":ee.ImageCollection('LANDSAT/LT04/C01/T1_SR'), 
			"L5": ee.ImageCollection('LANDSAT/LT05/C01/T1_SR'),
			"L7":ee.ImageCollection('LANDSAT/LE07/C01/T1_SR'),
			"L8":ee.ImageCollection("LANDSAT/LC08/C01/T1_SR")}
	
	if year < 1982 : 
		print("No landsat available")
	if year < 1993 and year > 1982:
		landsat = landsats['L4']
	if year < 2012 and year >1993:
		landsat = landsats['L5']
	if year < 2013 and year >=2012:
		landsat = landsats['L7']
	elif year >=2013: 
		landsat = landsats['L8']
	
	return landsat

def get_QA_bits(image, start, end, field_name):
	
	'''
	retrieve quality bits from landsat
	'''
	
	pattern = 0
	for i in range(start,end+1):
		pattern += 2**i
	return image.select([0], [field_name]).bitwiseAnd(pattern).rightShift(start)

def mask_quality(image):
	
	'''
	mask clouds and shoaws from landsat
	'''
	
	QA = image.select('pixel_qa')
	# Get the internal_cloud_algorithm_flag bit.
	shad = get_QA_bits(QA,3,3,'cloud_shadow')
	cloud = get_QA_bits(QA,5,5,'cloud')
	cirrus_detected = get_QA_bits(QA,9,9,'cirrus_detected')
	#Return an image masking out cloudy areas.
	return image.updateMask(shad.eq(0)).updateMask(cloud.eq(0).updateMask(cirrus_detected.eq(0))).unmask()

def get_data(year, start_doy, bounds):

	col1 = get_landsat(year).filter(ee.Filter.calendarRange(start_doy-8,start_doy+8, 'day_of_year')).filter(ee.Filter.calendarRange(year,year,'year'))
	col2 = col1.filterBounds(bounds).map(mask_quality)

	# landsat sr bandnames 
	bandnames = ["B1","B2", "B3", "B4", "B5", "B6", "B7"]

	ims = []
	for band in bandnames:
		ims.append(array_from_col(col2, band, 500, bounds))

	out = [np.nan_to_num(x) for x in ims]

	return out 


# main routine
def main():

	# Set the Study area (upoad kml to google fusion table, use the DocID in the ft string below)
	padded_cn = str(county_number).zfill(3)
	area = (ee.FeatureCollection('ft:1QPasan0i6O9uUlcYkjqj91D7mbnhTZCmzS4t7t_g').filter(ee.Filter().eq('id', padded_cn)))
	bounds = area.geometry().bounds()

	# Set the study years
	years = [x for x in range(2001, 2017)]
	band_names = ["B1","B2", "B3", "B4", "B5", "B6", "B7"]
	
	for year in years:
		print("Processing: " + str(year))

		# Setup outfiles like: /landsat/county_year_doy_band.txt
		outdir = os.path.join(os.getcwd(),"landsat")

		if not os.path.exists(outdir):
			os.mkdir(outdir)

		filepath = os.path.join(outdir,"_".join([str(county_number),str(year), str(start_doy)]))
		outpaths = [filepath+x+".txt" for x in band_names]

		# call earth engine to get the data 
		outdata = get_data(year, start_doy, bounds)
		
		# save
		for i in range(len(outdata)):
			np.savetxt(outpaths[i], outdata[i])

		print(" === Done ===")
	
main()
# CVWBM main readme

<img src="https://latex.codecogs.com/svg.latex?\Delta S_{gw} = P - ET - Q_{s} - Q_{g} \Delta S_{sm}" title="\Delta S = P - ET - Q_{s} - Q_{g}" />

![](https://ca.water.usgs.gov/projects/central-valley/images/ca3449_cover1.png)

### Lots of materials about the project
1. Data sources table
2. Time Periods Table
3. Links to all the data
4. etc...


# Satellite Data

## ET 

<sub>
    
|Product   			| ImageCollection 									| Var Name  	| Scaling Factor   	|   Availability	|   URL		|
|---				|---												|---			|---				|---				|---		|		
| modis   			| ee.ImageCollection('MODIS/006/MOD16A2') 			| "ET"  		|0.1   				|   2001 - Present	| [link](https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MOD16A2)  	|
| terraclimate  	| ee.ImageCollection('IDAHO_EPSCOR/TERRACLIMATE')  	| "aet"  		|0.1  				|   1980 - Present	| [link](https://developers.google.com/earth-engine/datasets/catalog/IDAHO_EPSCOR_GRIDMET)		|	
| gridmet eto 		| ee.ImageCollection('IDAHO_EPSCOR/GRIDMET')	  	| "eto"		  	|1				   	|   1980 - Present	| [link](https://developers.google.com/earth-engine/datasets/catalog/IDAHO_EPSCOR_GRIDMET)		|
| gridmet etr 		| ee.ImageCollection('IDAHO_EPSCOR/GRIDMET')	  	| "etr"		  	|1				   	|   1980 - Present	| [link]()		|
| fldas 	| ee.ImageCollection('NASA/FLDAS/NOAH01/C/GL/M/V001')   	| "Evap_tavg"  	|86400  			|   1980 - Present	| [link](https://developers.google.com/earth-engine/datasets/catalog/NASA_FLDAS_NOAH01_C_GL_M_V001)		|
| nldas  	| ee.ImageCollection('NASA/NLDAS/FORA0125_H002')  	|'potential_evaporation'|1		   			|   1980 - Present	| [link]		|
| gldas  	| ee.ImageCollection('NASA/GLDAS/V021/NOAH/G025/T3H' 		| "Evap_tavg"  	|86400  		  	|   1980 - Present	| [link]()		

</sub>

## P 

<sub>

|Product   			| ImageCollection 									| Var Name  		| Scaling Factor   	|   Availability	|   URL		|
|---				|---												|---				|---				|---				|--- 		|
| trmm   			| ee.ImageCollection('TRMM/3B43V7')					| "precipitation"  	|720  				| 2000 - Present  	|			|
| prism  			| ee.ImageCollection('OREGONSTATE/PRISM/AN81m')  	| "ppt"  			|1  				| 1980 - Present 	|			|
| chirps 			| ee.ImageCollection('UCSB-CHG/CHIRPS/PENTAD')	  	| "precipitation"	|1				   	|   				|			|
| persiann			| ee.ImageCollection('NOAA/PERSIANN-CDR')	  		| "precipitation"	|1				   	|   				|			|
| daymet 			| ee.ImageCollection('NASA/ORNL/DAYMET_V3')   		| "prcp"  			|1  				|   				|			|

</sub>

## SM 

<sub>
    
Product   			| ImageCollection 									| Var Name  		| Scaling Factor   	|   Availability	| URL		|
|---				|---												|---				|---				|---				|---
| terraclimate   	| ee.ImageCollection('IDAHO_EPSCOR/TERRACLIMATE')  	| "soil"  			|0.1  				|   1980 - Present	| []		|
</sub>

</sub>

```

asdf
###################
##### ET data #####
###################

# https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MOD16A2
data['modis'] = [ee.ImageCollection('MODIS/006/MOD16A2'), "ET", 0.1]

# https://developers.google.com/earth-engine/datasets/catalog/NASA_GLDAS_V020_NOAH_G025_T3H
data['gldas'] = [ee.ImageCollection('NASA/GLDAS/V021/NOAH/G025/T3H'), 'Evap_tavg', 86400]

# https://developers.google.com/earth-engine/datasets/catalog/NASA_NLDAS_FORA0125_H002
data['nldas'] = [ee.ImageCollection('NASA/NLDAS/FORA0125_H002'), 'potential_evaporation', 1]

# https://developers.google.com/earth-engine/datasets/catalog/IDAHO_EPSCOR_TERRACLIMATE
data['terra'] = [ee.ImageCollection('IDAHO_EPSCOR/TERRACLIMATE'), "aet", 0.1]

# https://developers.google.com/earth-engine/datasets/catalog/IDAHO_EPSCOR_GRIDMET
data['gmet_etr'] = [ee.ImageCollection('IDAHO_EPSCOR/GRIDMET'), "etr", 1]
data['gmet_eto'] = [ee.ImageCollection('IDAHO_EPSCOR/GRIDMET'), "eto", 1]

# https://developers.google.com/earth-engine/datasets/catalog/NASA_FLDAS_NOAH01_C_GL_M_V001
data['fldas'] = [ee.ImageCollection('NASA/FLDAS/NOAH01/C/GL/M/V001'), "Evap_tavg", 86400*24 ]


###################
##### P data ######
###################

data['trmm']  =  [ee.ImageCollection('TRMM/3B43V7'), "precipitation", 720]
data['prism'] = [ee.ImageCollection("OREGONSTATE/PRISM/AN81m"), "ppt", 1]
data['chirps'] = [ee.ImageCollection('UCSB-CHG/CHIRPS/PENTAD'), "precipitation", 1]
data['persia'] = [ee.ImageCollection("NOAA/PERSIANN-CDR"), "precipitation", 1]

https://developers.google.com/earth-engine/datasets/catalog/NASA_ORNL_DAYMET_V3
data['dmet'] = [ee.ImageCollection('NASA/ORNL/DAYMET_V3'), "prcp", 1]

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

data['sm1'] = [ee.ImageCollection("NASA/FLDAS/NOAH01/C/GL/M/V001"), "SoilMoi00_10cm_tavg", 86400*24]
data['sm2'] = [ee.ImageCollection("NASA/FLDAS/NOAH01/C/GL/M/V001"), "SoilMoi10_40cm_tavg", 86400*24]
data['sm3'] = [ee.ImageCollection("NASA/FLDAS/NOAH01/C/GL/M/V001"), "SoilMoi40_100cm_tavg", 86400*24]
data['sm4'] = [ee.ImageCollection("NASA/FLDAS/NOAH01/C/GL/M/V001"), "SoilMoi100_200cm_tavg", 86400*24]
```

# NHD Data for bounding boxes and watersheds: 
[https://www.usgs.gov/core-science-systems/ngp/national-hydrography/nhdplus-high-resolution](https://www.usgs.gov/core-science-systems/ngp/national-hydrography/nhdplus-high-resolution)
```
mkdir nhd
cd nhd

curl -o 1802.zip https://prd-tnm.s3.amazonaws.com/StagedProducts/Hydrography/NHD/HU4/HighResolution/Shape/NHD_H_1802_HU4_Shape.zip && mkdir 1802 && tar -xvf 1802.zip && mv Shape/ 1802/ && rm 1802.zip

curl -o 1803.zip https://prd-tnm.s3.amazonaws.com/StagedProducts/Hydrography/NHD/HU4/HighResolution/Shape/NHD_H_1803_HU4_Shape.zip && mkdir 1803 && tar -xvf 1803.zip && mv Shape/ 1803/ && rm 1803.zip

curl -o 1804.zip https://prd-tnm.s3.amazonaws.com/StagedProducts/Hydrography/NHD/HU4/HighResolution/Shape/NHD_H_1804_HU4_Shape.zip && mkdir 1804 && tar -xvf 1804.zip && mv Shape/ 1804/ && rm 1804.zip

```

# CDEC Reservoirs Data: 
[http://cdec.water.ca.gov/misc/monthlyStations.html](http://cdec.water.ca.gov/misc/monthlyStations.html) 


# CVWBM main readme

<a href="https://www.codecogs.com/eqnedit.php?latex=\Delta&space;S&space;=&space;P&space;-&space;ET&space;-&space;Q_{s}&space;-&space;Q_{g}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\Delta&space;S&space;=&space;P&space;-&space;ET&space;-&space;Q_{s}&space;-&space;Q_{g}" title="\Delta S = P - ET - Q_{s} - Q_{g}" /></a>

![](https://ca.water.usgs.gov/projects/central-valley/images/ca3449_cover1.png)

# Satellite Data

## ET 

<sub>
    
|Product   			| ImageCollection 									| Var Name  	| Scaling Factor   	|   Availability	|   URL		|
|---				|---												|---			|---				|---				|---		|		
| modis   			| ee.ImageCollection('MODIS/006/MOD16A2') 			| "ET"  		|0.1   				|   2001 - Present	| [link](https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MOD16A2)  	|
| terraclimate  	| ee.ImageCollection('IDAHO_EPSCOR/TERRACLIMATE')  	| "aet"  		|0.1  				|   1958 - Present	| [link](https://developers.google.com/earth-engine/datasets/catalog/IDAHO_EPSCOR_TERRACLIMATE)		|	
| gridmet eto 		| ee.ImageCollection('IDAHO_EPSCOR/GRIDMET')	  	| "eto"		  	|1				   	|   1980 - Present	| [link](https://developers.google.com/earth-engine/datasets/catalog/IDAHO_EPSCOR_GRIDMET)		|
| gridmet etr 		| ee.ImageCollection('IDAHO_EPSCOR/GRIDMET')	  	| "etr"		  	|1				   	|   1980 - Present	| [link](https://developers.google.com/earth-engine/datasets/catalog/IDAHO_EPSCOR_GRIDMET)		|
| fldas 	| ee.ImageCollection('NASA/FLDAS/NOAH01/C/GL/M/V001')   	| "Evap_tavg"  	|86400  			|   1980 - Present	| [link](https://developers.google.com/earth-engine/datasets/catalog/NASA_FLDAS_NOAH01_C_GL_M_V001)		|
| nldas  	| ee.ImageCollection('NASA/NLDAS/FORA0125_H002')  	|'potential_evaporation'|1		   			|   1980 - Present	| [link](https://developers.google.com/earth-engine/datasets/catalog/NASA_NLDAS_FORA0125_H002)		|
| gldas  	| ee.ImageCollection('NASA/GLDAS/V021/NOAH/G025/T3H' 		| "Evap_tavg"  	|86400  		  	|   1980 - Present	| [link](https://developers.google.com/earth-engine/datasets/catalog/NASA_GLDAS_V020_NOAH_G025_T3H)		

</sub>

## P 

<sub>

|Product   			| ImageCollection 									| Var Name  		| Scaling Factor   	|   Availability	|   URL		|
|---				|---												|---				|---				|---				|--- 		|
| trmm   			| ee.ImageCollection('TRMM/3B43V7')					| "precipitation"  	|720  				| 1998 - Present  	|[link](https://developers.google.com/earth-engine/datasets/catalog/TRMM_3B43V7)			|
| prism  			| ee.ImageCollection('OREGONSTATE/PRISM/AN81m')  	| "ppt"  			|1  				| 1895 - Present 	|[link](https://developers.google.com/earth-engine/datasets/catalog/OREGONSTATE_PRISM_AN81m)			|
| chirps 			| ee.ImageCollection('UCSB-CHG/CHIRPS/PENTAD')	  	| "precipitation"	|1				   	| 1981 - Present 	|[link](https://developers.google.com/earth-engine/datasets/catalog/UCSB-CHG_CHIRPS_PENTAD)			|
| persiann			| ee.ImageCollection('NOAA/PERSIANN-CDR')	  		| "precipitation"	|1				   	| 1983 - Present 	|[link](https://developers.google.com/earth-engine/datasets/catalog/NOAA_PERSIANN-CDR)			|
| daymet 			| ee.ImageCollection('NASA/ORNL/DAYMET_V3')   		| "prcp"  			|1  				| 1980 - Present 	|[link](https://developers.google.com/earth-engine/datasets/catalog/NASA_ORNL_DAYMET_V3)	|

</sub>

## SM 

<sub>
    
Product   			| ImageCollection 									| Var Name  		| Scaling Factor   	|   Availability	| URL		|
|---				|---												|---				|---				|---				|---
| terraclimate   	| ee.ImageCollection('IDAHO_EPSCOR/TERRACLIMATE')  	| "soil"  			|0.1  				|   1958 - Present	| [link](https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MOD16A2)		|
</sub>

</sub>

## Implementation
```
These data are implimented in code like this: 

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


###################
##### P data ######
###################

data['trmm']  =  [ee.ImageCollection('TRMM/3B43V7'), "precipitation", 720]
data['prism'] = [ee.ImageCollection("OREGONSTATE/PRISM/AN81m"), "ppt", 1]
data['chirps'] = [ee.ImageCollection('UCSB-CHG/CHIRPS/PENTAD'), "precipitation", 1]
data['persia'] = [ee.ImageCollection("NOAA/PERSIANN-CDR"), "precipitation", 1]
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

# Watershed Data - National Hydrographic Dataset 
[https://www.usgs.gov/core-science-systems/ngp/national-hydrography/nhdplus-high-resolution](https://www.usgs.gov/core-science-systems/ngp/national-hydrography/nhdplus-high-resolution)
```
mkdir nhd
cd nhd

curl -o 1802.zip https://prd-tnm.s3.amazonaws.com/StagedProducts/Hydrography/NHD/HU4/HighResolution/Shape/NHD_H_1802_HU4_Shape.zip && mkdir 1802 && tar -xvf 1802.zip && mv Shape/ 1802/ && rm 1802.zip

curl -o 1803.zip https://prd-tnm.s3.amazonaws.com/StagedProducts/Hydrography/NHD/HU4/HighResolution/Shape/NHD_H_1803_HU4_Shape.zip && mkdir 1803 && tar -xvf 1803.zip && mv Shape/ 1803/ && rm 1803.zip

curl -o 1804.zip https://prd-tnm.s3.amazonaws.com/StagedProducts/Hydrography/NHD/HU4/HighResolution/Shape/NHD_H_1804_HU4_Shape.zip && mkdir 1804 && tar -xvf 1804.zip && mv Shape/ 1804/ && rm 1804.zip

```

# Reservoirs Data - California Data Exchange Center
[http://cdec.water.ca.gov/misc/monthlyStations.html](http://cdec.water.ca.gov/misc/monthlyStations.html) 


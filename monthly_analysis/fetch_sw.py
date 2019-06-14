import os
import shapely
import io
import requests
import datetime
import urllib.request

import geopandas as gp
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
    
from shapely.geometry import Point
from shapely.ops import cascaded_union, unary_union
from shapely.geometry import Point

from climata.usgs import DailyValueIO

def get_streamflow(huc8):
    
    '''
    call climata API supplying huc8 argument to get each gaging station within each basin 
    '''
    
    data =  DailyValueIO (
            start_date="2001-01-01", 
            end_date="2018-01-01",
            basin=huc8,
            parameter="00060",
            )
    
    qs = []
    ds = []
    lats = []
    lons = []
    ids = []

    for series in data:
        values = []
        dates = []
        lats.append(series.latitude)
        lons.append(series.longitude)
        ids.append(series.site_code)

        for row in series.data:
            values.append(row.value)
            dates.append(row.date)

        qs.append(values)
        ds.append(dates)
    
    geometry = [Point(xy) for xy in zip(lons, lats)]
    df = pd.DataFrame(geometry)
    crs = {'init': 'epsg:4326'}
    gdf = gp.GeoDataFrame(df, crs=crs, geometry=geometry)
    
    return gdf, qs, ds, ids

def find_hucs(dir, huc_str):
    shp_files = [os.path.join(dir,x) for x in os.listdir(dir) if huc_str in x if "xml" not in x]
    return shp_files

def main_streamflow():
    # find the Huc 8s that climata accepts to query. 
    huc_order = "8"
    huc_str = "WBDHU{}.shp".format(huc_order)

    huc4s = ["1802", "1803", "1804"]
    hu4_dirs = [os.path.join("../nhd", x, "Shape") for x in os.listdir("../nhd") if "." not in x]


    gdfs = []

    for i in hu4_dirs:
        gdfs.append(gp.read_file(find_hucs(i, huc_str)[0]))

    hu4 = gp.GeoDataFrame(pd.concat(gdfs, ignore_index=True))

    gdfs = []
    qs = []
    ds = []
    ids = []

    for i in hu4['HUC8']:
        print ("processing " + i)
        gdf, q, d, i = get_streamflow(i)
        gdfs.append(gdf)
        qs.append(q)
        ds.append(d)
        ids.append(i)

    # Un-nest the lists

    q_f = np.array([item for sublist in qs for item in sublist])
    d_f = np.array([item for sublist in ds for item in sublist])
    ids_f = [item for sublist in ids for item in sublist]

    # Stations from Xiao:
    stations = [11446500, 11376550, 11423800, 11384000, 11390000 ,11451760,
                11372000, 11335000, 11376000, 11374000, 11383500, 11329500,
                11211300, 11424500, 11379500, 11407150, 11257500, 11209900,
                11192950, 11251600, 11225000, 11270900, 11381500, 11221700,
                11325500, 11384350, 11454000, 11370500, 11251000, 11302000,
                11388000, 11382000, 11289650, 11199500, 11421000]

    # Stations I added in last row 
    stations = [11446500, 11376550, 11423800, 11384000, 11390000 ,11451760,
                11372000, 11335000, 11376000, 11374000, 11383500, 11329500,
                11211300, 11424500, 11379500, 11407150, 11257500, 11209900,
                11192950, 11251600, 11225000, 11270900, 11381500, 11221700,
                11325500, 11384350, 11454000, 11370500, 11251000, 11302000, 
                11388000, 11382000, 11289650, 11199500, 11421000, 
               
                11208818, 11204100, 11200800, 11218400, 11289000, 11323500
               ]

    # The CA Aqueduct takes water out of the CV: 
    stations_out = ["11109396"]

    stations = [str(x) for x in stations]

    # Make a gdf of the stations and join the lists of q, ID and dates
    stations_gdf = gp.GeoDataFrame(pd.concat(gdfs, ignore_index=True, sort = False))
    stations_gdf['Q'] = q_f
    stations_gdf['ID'] = ids_f
    stations_gdf['dates'] = d_f

    # Filter the master GDF for the relevant stations 
    inflow = stations_gdf[stations_gdf['ID'].isin(stations)]
    outflow = stations_gdf[stations_gdf.ID == "11109396"]

    q_arrs = []
    d_arrs = []

    for i in range(len(inflow)):
        q_arrs.append(inflow.iloc[i].Q)
        d_arrs.append(inflow.iloc[i].dates)

    dfs = []

    for i,x in enumerate(q_arrs):
        df = pd.DataFrame({'q':q_arrs[i],'date':d_arrs[i]})
        df = df.set_index(pd.DatetimeIndex(df['date']))
        df.drop(['date'], axis = 1, inplace = True)
        monthly =  df.resample('M').sum()
        dfs.append(monthly)

    # The outflow here is the CA Aqueduct. If running for another area besides CV, need to add outflows here: 

    # If there are more outflows, need to make this a loop 
    q_out = outflow.iloc[0].Q
    q_out_d = outflow.iloc[0].dates

    df_o = pd.DataFrame({'q_out':q_out,'date':q_out_d})
    df_o = df_o.set_index(pd.DatetimeIndex(df_o['date']))
    df_o.drop(['date'], axis = 1, inplace = True)
    monthly_out =  df_o.resample('M').sum() * -1

    dfs.append(monthly_out)

    t = pd.concat(dfs, axis=1)
    sums = t.sum(axis=1, skipna=True)

    # cfs to km^3/month
    sums = sums * 0.0283168 * 1e-9 * 86400
    sums.to_csv("Qs_monthly.csv")
    print("Fetched Discharge Data")

def main_reservoirs():
    # Load cdec reservoirs
    data_dir = "../data"
    files = [os.path.join(data_dir, x) for x in os.listdir(data_dir) if "cdec" in x]
    df = pd.read_csv(files[0])

    # Make a gdf
    lats = df.Latitude
    longs = df.Longitude
    crs = {"init":"epsg:4326"}
    geom = [Point(xy) for xy in zip(longs,lats)]
    gdf = gp.GeoDataFrame(df, crs = crs, geometry = geom)

    # Read the HU4s and merge
    s1 = gp.read_file("../shape/1804_4.shp")
    s2 = gp.read_file("../shape/1803_4.shp")
    s3 = gp.read_file("../shape/1802_4.shp")

    ssjt = gp.GeoDataFrame(pd.concat([s1,s2,s3], ignore_index=True))
    ssjt.crs = crs

    # Read cv
    cv = gp.read_file("../shape/cv.shp")

    # Make small buffer around cv
    buffered_cv = cv.buffer(0.005)
    buf = gp.GeoDataFrame(gp.GeoSeries(buffered_cv))
    buf = buf.rename(columns={0:'geometry'}).set_geometry('geometry')
    buf.crs = crs

    # Filter
    within_cv = gp.sjoin(gdf, buf, how='inner', op='intersects')
    within_ssjt = gp.sjoin(gdf, ssjt, how='inner', op='intersects')

    # Download Reservoir Storage (SensorNums = 15) data by query str:
    start = datetime.datetime(2001, 1, 1)
    end = datetime.datetime(2019, 1, 1)
    dt_idx = pd.date_range(start,end, freq='M')

    ssjt_data = {}

    for i in within_ssjt.ID:
        print("processing " + i )
        url = "https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations={}&SensorNums=15&dur_code=M&Start=2001-01-01&End=2018-12-01".format(i)
        urlData = requests.get(url).content
        df = pd.read_csv(io.StringIO(urlData.decode('utf-8')))
        
        if df.empty:
            pass
        else:
            ssjt_data[i] = df

    cv_data={}

    for i in within_cv.ID:
        print("processing " + i )
        url = "https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations={}&SensorNums=15&dur_code=M&Start=2001-01-01&End=2018-12-01".format(i)
        urlData = requests.get(url).content
        df = pd.read_csv(io.StringIO(urlData.decode('utf-8')))
        
        if df.empty:
            pass
        else:
            cv_data[i] = df

    ssjt_s = []
    ssjt_df = {}

    for k,v in ssjt_data.items():
        ssjt_s.append(pd.to_numeric(data[k].VALUE, errors = "coerce"))
        ssjt_df[k] = ssjt_data[k].VALUE

    cv_s = []
    cv_df = {}
    for k,v in cv_data.items():
        cv_s.append(pd.to_numeric(data[k].VALUE, errors = "coerce"))
        ssjt_df[k] = data[k].VALUE

    ssjt_ts = np.nansum(np.column_stack(ssjt_s), axis = 1)
    cv_ts = np.nansum(np.column_stack(cv_s), axis = 1)

    out = pd.DataFrame([ssjt_ts*1.23348e-6, cv_ts*1.23348e-6]).T
    out.columns = ['ssjt_res_s_km3',"cv_res_s_km3"]
    out.index = dt_idx
    out.head()

    out.to_csv("reservoir_storage.csv")
    print("Fetched Reservoir Storage Data") 

def main():
    # main_streamflow()
    main_reservoirs()

if __name__ == "__main__":
    main()
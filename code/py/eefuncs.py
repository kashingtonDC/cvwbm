import numpy as np
import pandas as pd
import ee

# Earth Engine Processing functions

def df_from_imcol(imcol):

    '''
    Takes an earth engine class (ee.ImageCollection) and converts it to a pandas df
    (sanity check)
    '''

    df = pd.DataFrame(imcol, columns = imcol[0])
    df = df[1:]

    return(df)

def array_from_df(df, col_name):

    '''
    Takes a pandas data frame and converts to a numpy array which can be plotted, ingested, etc
    '''
    
    # get data from df as arrays
    lons = np.array(df.longitude)
    lats = np.array(df.latitude)
    data = np.array(df[col_name]) # Name of column
    
    # get the unique coordinates
    uniqueLats = np.unique(lats)
    uniqueLons = np.unique(lons)

    # get number of columns and rows from coordinates
    ncols = len(uniqueLons)    
    nrows = len(uniqueLats)

    # determine pixelsizes
    ys = uniqueLats[1] - uniqueLats[0] 
    xs = uniqueLons[1] - uniqueLons[0]

    # create an array with dimensions of image
    arr = np.zeros([nrows, ncols], np.float32)

    # fill the array with values
    counter =0
    for y in range(0,len(arr),1):
        for x in range(0,len(arr[0]),1):
            if lats[counter] == uniqueLats[y] and lons[counter] == uniqueLons[x] and counter < len(lats)-1:
                counter+=1
                arr[len(uniqueLats)-1-y,x] = data[counter] # start from lower left corner
    
    return arr

def array_from_col(col,band,res,bounds):
    
    '''
    Transform an ee.ImageCollection class to a numpy array
    
    Args: 
    
    col: ee.ImageColletion ex 'LANDSAT/LT04/C01/T1_SR'
    band: string, ex "B1" '
    res: int, ex: 30
    

    #     start = ee.Date.fromYMD(year,month,day)
    #     end = start.advance(1,'month')
    #     col = get_landsat(year).filterBounds(area).filterDate(start, end).map(mask_quality)
        

    '''

    # get the lat lon and add the band and scale by the appropriate factor (0.0001 for landsat)
    band_name = col.select(band).median()
    latlon = ee.Image.pixelLonLat().addBands(band_name).multiply(0.0001)

    # apply reducer to list
    latlon = latlon.reduceRegion(
      reducer=ee.Reducer.toList(),
      geometry=bounds,
      maxPixels=1e13,
      scale=res)
    
    data = np.array((ee.Array(latlon.get(band)).getInfo()))
    lats = np.array((ee.Array(latlon.get("latitude")).getInfo()))
    lons = np.array((ee.Array(latlon.get("longitude")).getInfo()))
    
    arr = array_from_coords(data,lats,lons)
    
    return (arr)

def array_from_coords(data,lats,lons):
    
    '''
    Return a numpy array (ie cartesian product) from lats, lons, and data values
    '''
    
    # get the unique coordinates
    uniqueLats = np.unique(lats)
    uniqueLons = np.unique(lons)

    # get number of columns and rows from coordinates
    ncols = len(uniqueLons)    
    nrows = len(uniqueLats)

    # determine pixelsizes
    ys = uniqueLats[1] - uniqueLats[0] 
    xs = uniqueLons[1] - uniqueLons[0]

    # create an array with dimensions of image
    arr = np.zeros([nrows, ncols], np.float32) #-9999

    # fill the array with values
    counter =0
    for y in range(0,len(arr),1):
        for x in range(0,len(arr[0]),1):
            if lats[counter] == uniqueLats[y] and lons[counter] == uniqueLons[x] and counter < len(lats)-1:
                counter+=1
                arr[len(uniqueLats)-1-y,x] = data[counter] 
                
    return arr

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

def gen_polys(geometry, dx, dy):
    '''
    Split up an area in to a lot of polygons so you can query earth engine. 
    '''
    bounds = ee.Geometry(geometry).bounds()
    coords = ee.List(bounds.coordinates().get(0))
    ll = ee.List(coords.get(0))
    ur = ee.List(coords.get(2))
    xmin = ll.get(0)
    xmax = ur.get(0)
    ymin = ll.get(1)
    ymax = ur.get(1)

    latlist = ee.List.sequence(ymin, ymax, dx)
    lonlist = ee.List.sequence(xmin, xmax, dy)
    
    for lat in latlist.getInfo():
        def make_rect(lat):
            lattemp = ee.Number(lat)
            for lon in lonlist.getInfo():
                lontemp = ee.Number(lon)
                uplattemp = lattemp.subtract(dy)
                lowlontemp = lontemp.add(dx)
                return ee.Feature(ee.Geometry.Polygon([[lontemp, lattemp],
                                                        [lowlontemp, lattemp],
                                                        [lowlontemp, uplattemp],
                                                        [lontemp, uplattemp]]))
    s = latlist.map(make_rect).flatten()       
    
    return s 


'''
Stats functions - do these belong elsewhere? 

'''

def moving_window(image, window_w = 8, window_h = 8):
    '''
    Creates a list of 8x8 chips of an image
    '''
    w, h = image.shape
    w_new, h_new = w - window_w + 1, h - window_h + 1
    scans = []
    for i in range(w_new):
        for j in range(h_new):
            scans.append(image[i:i+window_w, j:j+window_h])
    return scans

def get_pc(image):
    
    '''
    apply an 8x8 pixel (2x2 km) moving window to modis imagery, 
    return the top3 principal components as images and corresponding normalized eigenvectors (explaining % of variance)

    '''
    
    window_h = 8
    window_w = 8
    h,w = image.shape
    nrows,ncols = h+1-window_h, w+1-window_w
    mw = moving_window(image)
    
    flat = [np.array(x).flatten() for x in mw] # Flatten the 5x5 windows 
    array = np.vstack(flat) # stack to a big array

    pca = decomposition.PCA()
    tf = pca.fit_transform(array)
    
    # return first 3 pcs and first 3 eigenvalues
    pc1 = tf.T[0].reshape(nrows,ncols)
    ev1 = pca.fit(array).explained_variance_ratio_[0]
    pc2 = tf.T[1].reshape(nrows,ncols)
    ev2 = pca.fit(array).explained_variance_ratio_[1]
    pc3 = tf.T[2].reshape(nrows,ncols)
    ev3 = pca.fit(array).explained_variance_ratio_[2]

    return ([[pc1, pc2, pc3] , [ev1, ev2, ev3]])
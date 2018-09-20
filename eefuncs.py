
'''
This module contains functions to handle conversion between ee.ImageCollection class and numpy.ndarray class

Processing flow is:
ee.ImageCollection --> pandas.DataFrame --> numpy.ndarray

'''

# Libs

import numpy as np
import pandas as pd

# Functions

def df_from_imcol(imcol):

    '''
    Takes an earth engine class (ee.ImageCollection) and converts it to a pandas df
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
    data = np.array(df[col_name])
    #data = np.array([np.float(x) for x in np.array(df[col_name])]) # Ensure that the data are converted to float
    
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
    arr = np.zeros([nrows, ncols], np.float16)

    # fill the array with values
    counter =0
    for y in range(0,len(arr),1):
        for x in range(0,len(arr[0]),1):
            if lats[counter] == uniqueLats[y] and lons[counter] == uniqueLons[x] and counter < len(lats)-1:
                counter+=1
                arr[len(uniqueLats)-1-y,x] = data[counter] # start from lower left corner

    arr = arr.astype(np.float)
    
    return (arr)
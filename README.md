### Geospatial Build

Aakash Ahamed
Stanford University
September, 2018

### This readme contains instructions to build a geospatial data processing environment which includes:

- Python 3.6
- Earth Engine
- GDAL, Geopandas, Rasterio, Fiona, other geospatial libs
- Numpy, scipy, pandas, other useful python libs
- Tensorflow, Keras, sklearn, other ML libs

### Download Anaconda3

### Clone this repository
`git clone **this repo**`

### Create a new conda virtual environment from the env.yml file
`conda env create -f env.yml`
 
This will take a while to install the packages. The name of the environment is 'gis', which is specified in the first line of the env.yml file

Activate the environment:
`source activate gis`


Then (gis) will be prepended to your path. Note that for newer versions of conda, this command may be `conda activate gis`

Next run:
`pip install google-api-python-client`
`conda install cython`
`pip install pyCrypto`
`pip install earthengine-api`
`pip install oauth2client`

Now install Tensowflow backend:
`pip install --ignore-installed --upgrade https://storage.googleapis.com/tensorflow/mac/cpu/tensorflow-1.10.1-py3-none-any.whl`

And Keras:
`pip install keras`

Test everything:
```
python
import ee
import gdal
import tensorflow as tf
import keras
ee.Initialize() # Supply your EE credentials
```
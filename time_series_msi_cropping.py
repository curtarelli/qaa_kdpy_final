'''
@author: Victor Pedroso Curtarelli
-------------------------------------------------------------------------------


'''
import os
import gc
import glob
from osgeo import gdal
from gdalconst import *
import pandas as pd
import geopandas as gpd
import numpy as np
import numpy.ma as ma
import math as mt
import datetime as dt
import scipy.ndimage

from shapely.geometry import mapping
import rasterio as rio
from rasterio.plot import plotting_extent
from rasterio.mask import mask
import earthpy as et
import earthpy.spatial as es
import earthpy.plot as ep

import seaborn as sns
from matplotlib import pyplot as plt

# Prettier plotting with seaborn
sns.set(font_scale = 1.5)

##  Diretório onde se encontram as imagens
dir_i = r'C:\Imagens'
dir_shp = r'C:\Imagens\study_area\study_area_utm.shp'

crop_extent = gpd.read_file(dir_shp)
print('crop extent crs: ', crop_extent.crs)

'''
# Plot the crop boundary layer
# Note this is just an example so you can see what it looks like
# You don't need to plot this layer in your homework!
fig, ax = plt.subplots(figsize=(6, 6))

crop_extent.plot(ax=ax)

ax.set_title("Shapefile Crop Extent",
             fontsize=16)
'''

os.chdir(dir_i)
list_i = glob.glob('20*')

# i = list_i[0]

for i in list_i:
    dir_ii = dir_i + '\\' + i
    os.chdir(dir_ii)
    print('Entrou em ' + i)
    
    list_ii = glob.glob('*')
    
    # ii = list_ii[0]
    
    for ii in list_ii:        
        dir_iii = dir_ii + '\\' + ii
        
        os.chdir(dir_iii)
        print('Entrou em ' + i + '\\' + ii)
        
        dir_crop = dir_iii + '\\crop'
        
        if not os.path.exists(dir_crop):
            os.mkdir(dir_crop)

        dir_crop = dir_iii + '\\crop'
        
        list_iii = glob.glob('*.tif')
        
        # iii = list_iii[1]
        
        for iii in list_iii:
            
            with rio.open(iii) as iii_tif:
                iii_tif_crop, iii_tif_crop_meta = es.crop_image(iii_tif, crop_extent)
            
            iii_tif_crop_affine = iii_tif_crop_meta['transform']
            
            # Create spatial plotting extent for the cropped layer
            iii_tif_extent = plotting_extent(iii_tif_crop[0], iii_tif_crop_affine)
            
            # Plot your data
            '''
            ep.plot_bands(iii_tif_crop[0],
                          extent = iii_tif_extent,
                          cmap = 'Greys',
                          title = "Cropped Raster Dataset",
                          scale = False)
            plt.show()
            '''
            
            # Update with the new cropped affine info and the new width and height
            iii_tif_crop_meta.update({'transform': iii_tif_crop_affine,
                                   'height': iii_tif_crop.shape[1],
                                   'width': iii_tif_crop.shape[2],
                                   'nodata': int(255)})
            iii_tif_crop_meta
            
            crop_name = dir_crop + '\\crop_' + iii
            
            with rio.open(crop_name, 'w', **iii_tif_crop_meta) as ff:
                ff.write(iii_tif_crop[0], 1)



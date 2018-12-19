# -*- coding: utf-8 -*-
"""
Created on 12 December 2018

@author: Heiko Balzter

"""

#############################################################################
# read all classified images in a directory and a shape file
#   and make jpeg quicklook maps at different scales
# written for Python 3.6.4
#############################################################################

import matplotlib.pyplot as plt
plt.switch_backend('agg') # solves QT5 problem
import os, sys
from osgeo import gdal, gdalnumeric, ogr, osr
from skimage import io

gdal.UseExceptions()
io.use_plugin('matplotlib')

#############################################################################
# OPTIONS
#############################################################################
copyright = '© University of Leicester, 2018. ' #text to be plotted on the map
wd = '/scratch/clcr/shared/heiko/marque_de_com/images/' # working directory on Linux HPC
shapedir = '/scratch/clcr/shared/heiko/aois/' # this is where the shapefile is
datadir = wd + 'L2/'  # directory of Sentinel L2A data files in .SAFE format
classdir = wd + 'class/'  # directory of classified images
mapdir = wd + 'maps/' # directory for L2A maps
classmapdir = wd + 'classmaps/'  # directory for classified maps
shapefile = shapedir + 'marque.shp' # shapefile of test area
bands = ['B04_10m','B03_10m','B02_10m'] #corresponds to 10 m resolution Sentinel-2 bands Red, Green, Blue for image display
rosepath = '/home/h/hb91/PycharmProjects/pyeo/pyeo/' # location of compassrose.jpg on HPC


#############################################################################
# MAIN
#############################################################################

# go to working directory
os.chdir(wd)

# make a 'classmaps' directory (if it does not exist yet) for output files
if not os.path.exists(classmapdir):
    print("Creating directory: ", classmapdir)
    os.mkdir(classmapdir)

n = pyeo.map_all_class_images(classdir, id="Overview", cols=None, figsizex=12, figsizey=12, zoom=1, xoffset=0, yoffset=0) # overview map
print("Made "+str(n)+" maps.")
n = pyeo.map_all_class_images(classdir, id="ZoomOut", cols=None, figsizex=12, figsizey=12, zoom=2, xoffset=0, yoffset=0) # zoom out
print("Made "+str(n)+" maps.")
n = pyeo.map_all_class_images(classdir, id="ZoomIn", cols=None, figsizex=12, figsizey=12, zoom=0.1, xoffset=0, yoffset=0) # zoom in
print("Made "+str(n)+" maps.")
n = pyeo.map_all_class_images(classdir, id="MoveLeft", cols=None, figsizex=12, figsizey=12, zoom=0.1, xoffset=0, yoffset=-2500) # move left
print("Made "+str(n)+" maps.")

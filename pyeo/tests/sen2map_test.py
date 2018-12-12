# -*- coding: utf-8 -*-
"""
Created on 6 December 2018

@author: Heiko Balzter

"""

#############################################################################
# read all Sentinel-2 scenes in a directory and a shape file
#   and make RGB quicklook maps at different scales
# written for Python 3.6.4
#############################################################################

# When you start the IPython Kernel, launch a graphical user interface (GUI) loop:
#   %matplotlib

from cartopy.io.shapereader import Reader
from cartopy.feature import ShapelyFeature, BORDERS
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import datetime
import json
import matplotlib.image as im
import matplotlib.lines as mlines
import matplotlib.patches as patches
import matplotlib.pyplot as plt
plt.switch_backend('agg') # solves QT5 problem
import numpy as np
import os, sys
from os import listdir
from os.path import isfile, isdir, join
from osgeo import gdal, gdalnumeric, ogr, osr
from skimage import io
import subprocess

gdal.UseExceptions()
io.use_plugin('matplotlib')

# suppress warnings when dividing by zero or nan
np.seterr(divide='ignore', invalid='ignore')

# suppress matplotlib depracation warning
#warnings.filterwarnings("ignore",category=matplotlib.cbook.mplDeprecation)

# The pyplot interface provides 4 commands that are useful for interactive control.
# plt.isinteractive() returns the interactive setting True|False
# plt.ion() turns interactive mode on
# plt.ioff() turns interactive mode off
# plt.draw() forces a figure redraw


#############################################################################
# OPTIONS
#############################################################################
copyright = '© University of Leicester, 2018. ' #text to be plotted on the map
wd = '/scratch/clcr/shared/heiko/marque_de_com/images/' # working directory on Linux HPC
shapedir = '/scratch/clcr/shared/heiko/aois/' # this is where the shapefile is
datadir = wd + 'L2/'  # directory of Sentinel L2A data files in .SAFE format
mapdir = wd + 'maps/'  # directory of Sentinel L1C data files in .SAFE format
shapefile = shapedir + 'marque.shp' # shapefile of test area
bands = ['B04_10m','B03_10m','B02_10m'] #corresponds to 10 m resolution Sentinel-2 bands Red, Green, Blue for image display
rosepath = '/home/h/hb91/PycharmProjects/pyeo/pyeo/' # location of compassrose.jpg on HPC

#############################################################################
# FUNCTION DECLARATIONS
#############################################################################

def blank_axes(ax):
    """
    blank_axes:  blank the extraneous spines and tick marks for an axes

    Input:
    ax:  a matplotlib Axes object

    Output: None
    """

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.yaxis.set_ticks_position('none')
    ax.xaxis.set_ticks_position('none')
    ax.tick_params(labelbottom='off', labeltop='off', labelleft='off', labelright='off', \
                   bottom='off', top='off', left='off', right='off')

# define functions to read/write floating point numbers from/to a text file
def read_floats(filename):
    with open(filename) as f:
        return [float(x) for x in f]
    f.close()

def write_floats(data, filename):
    file = open(filename, 'w')
    for item in data:
        file.write("%f\n" % item)
    file.close()

def get_gridlines(x0, x1, y0, y1, nticks):
    '''
    make neat gridline labels for map projections
        x0, x1 = minimum and maximum x positions in map projection coordinates
        y0, y1 = minimum and maximum y positions in map projection coordinates
        nticks = number of ticks / gridlines in x direction
        returns a numpy array with x and y tick positions
    '''
    # calculate length of axis
    lx = x1 - x0

    # count number of digits of axis lengths
    nlx = int(np.log10(lx) + 1)

    # divide lengths into segments and round to highest digit
    #   remove all but the highest digit
    ndigits = int(np.log10(lx / nticks))
    dx = int(lx / nticks / 10 ** ndigits)
    #   round to a single digit integer starting with 1, 2 or 5
    pretty = [1, 2, 5, 10] # pretty numbers for the gridlines
    d = [0, 0, 0, 0] # absolute differences between dx and pretty numbers
    d[:] = [abs(x - dx) for x in pretty]
    # find the index of the pretty number with the smallest difference to dx and then the number
    dx = pretty[np.argmin(d)]
    #   scale back up
    dx = dx * 10 ** ndigits
    # update number of digits in case pretty is 10
    ndigits = int(np.log10(dx))

    # find position of the first pretty gridline just outside the map area
    xs = int(x0 / 10 ** ndigits) * 10 ** ndigits

    # set x ticks positions
    xticks = np.arange(xs, x1 + dx -1, dx)
    #xticks = [x for x in xt if (x >= x0 and x <=x1)] # checks whether outside of map boundary, not needed

    # find position of the first pretty gridline just outside the map area
    ys = int(y0 / 10 ** ndigits) * 10 ** ndigits

    # set y ticks positions
    yticks = np.arange(ys, y1 + dx -1, dx)

    return xticks, yticks

# function to convert coordinates
def convertXY(xy_source, inproj, outproj):
    shape = xy_source[0, :, :].shape
    size = xy_source[0, :, :].size
    # the ct object takes and returns pairs of x,y, not 2d grids
    # so the the grid needs to be reshaped (flattened) and back.
    ct = osr.CoordinateTransformation(inproj, outproj)
    xy_target = np.array(ct.TransformPoints(xy_source.reshape(2, size).T))
    xx = xy_target[:, 0].reshape(shape)
    yy = xy_target[:, 1].reshape(shape)
    return xx, yy


# This function will convert the rasterized clipper shapefile to a mask for use within GDAL.
def imageToArray(i):
    """
    Converts a Python Imaging Library array to a
    gdalnumeric image.
    """
    a = gdalnumeric.fromstring(i.tostring(), 'b')
    a.shape = i.im.size[1], i.im.size[0]
    return a


def world2Pixel(geoMatrix, x, y):
    """
    Uses a gdal geomatrix (gdal.GetGeoTransform()) to calculate
    the pixel location of a geospatial coordinate
    """
    ulX = geoMatrix[0]
    ulY = geoMatrix[3]
    xDist = geoMatrix[1]
    yDist = geoMatrix[5]
    rtnX = geoMatrix[2]
    rtnY = geoMatrix[4]
    pixel = int((x - ulX) / xDist)
    line = int((ulY - y) / xDist)
    return (pixel, line)


def transformxy(s_srs, t_srs, xcoord, ycoord):
    """
    Transforms a point coordinate x,y from a source reference system (s_srs)
    to a target reference system (t_srs)
    """
    geom = ogr.Geometry(ogr.wkbPoint)
    geom.SetPoint_2D(0, xcoord, ycoord)
    geom.AssignSpatialReference(s_srs)
    geom.TransformTo(t_srs)
    return geom.GetPoint_2D()


def projectshape(inshp, outshp, t_srs):
    """
    Reprojects an ESRI shapefile from its source reference system
    to a target reference system (e.g. t_srs = 4326)
    filenames must include the full directory paths
    requires:
        from osgeo import ogr, osr
        import os
    """

    driver = ogr.GetDriverByName('ESRI Shapefile')  # get shapefile driver
    infile = driver.Open(inshp, 0)
    if infile is None:
        print('Could not open ' + inshp)
        sys.exit(1)  # exit with an error code
    inLayer = infile.GetLayer()  # get input layer
    inSpatialRef = inLayer.GetSpatialRef()  # get source spatial reference system
    # or input SpatialReference manually here
    #   inSpatialRef = osr.SpatialReference()
    #   inSpatialRef.ImportFromEPSG(2927)
    outSpatialRef = osr.SpatialReference()
    outSpatialRef.ImportFromEPSG(t_srs)
    # create the CoordinateTransformation
    coordTrans = osr.CoordinateTransformation(inSpatialRef, outSpatialRef)
    # create the output layer
    if os.path.exists(outshp):
        driver.DeleteDataSource(outshp)
    outDataSet = driver.CreateDataSource(outshp)
    outLayer = outDataSet.CreateLayer("basemap_" + str(t_srs), geom_type=ogr.wkbMultiPolygon)
    # add fields
    inLayerDefn = inLayer.GetLayerDefn()
    for i in range(0, inLayerDefn.GetFieldCount()):
        fieldDefn = inLayerDefn.GetFieldDefn(i)
        outLayer.CreateField(fieldDefn)
    # get the output layer's feature definition
    outLayerDefn = outLayer.GetLayerDefn()
    # loop through the input features
    inFeature = inLayer.GetNextFeature()
    while inFeature:
        # get the input geometry
        geom = inFeature.GetGeometryRef()
        # reproject the geometry
        geom.Transform(coordTrans)
        # create a new feature
        outFeature = ogr.Feature(outLayerDefn)
        # set the geometry and attribute
        outFeature.SetGeometry(geom)
        for i in range(0, outLayerDefn.GetFieldCount()):
            outFeature.SetField(outLayerDefn.GetFieldDefn(i).GetNameRef(), inFeature.GetField(i))
        # add the feature to the shapefile
        outLayer.CreateFeature(outFeature)
        # dereference the features and get the next input feature
        outFeature = None
        inFeature = inLayer.GetNextFeature()
    # Save and close the shapefiles
    inDataSet = None
    outDataSet = None
    # Try to open the output file to check it worked
    outfile = driver.Open(outshp, 0)
    if outfile is None:
        print('Failed to create ' + outshp)
        sys.exit(1)  # exit with an error code
    else:
        print('Reprojection of shapefile seems to have worked.')
    return None


def OpenArray(array, prototype_ds=None, xoff=0, yoff=0):
    #  this is basically an overloaded version of the gdal_array.OpenArray passing in xoff, yoff explicitly
    #  so we can pass these params off to CopyDatasetInfo
    ds = gdal.Open(gdalnumeric.GetArrayFilename(array))

    if ds is not None and prototype_ds is not None:
        if type(prototype_ds).__name__ == 'str':
            prototype_ds = gdal.Open(prototype_ds)
        if prototype_ds is not None:
            gdalnumeric.CopyDatasetInfo(prototype_ds, ds, xoff=xoff, yoff=yoff)
    return ds


def histogram(a, bins=range(0, 256)):
    """
    Histogram function for multi-dimensional array.
    a = array
    bins = range of numbers to match
    """
    fa = a.flat
    n = gdalnumeric.searchsorted(gdalnumeric.sort(fa), bins)
    n = gdalnumeric.concatenate([n, [len(fa)]])
    hist = n[1:] - n[:-1]
    return hist


def stretch(im, nbins=256, p=None, nozero=True):
    """
    Performs a histogram stretch on an ndarray image.
    im = image
    nbins = number of histogram bins
    p = percentile to be removed at the bottom and top end of the histogram (0-100)
    nozero = remove zero values from histogram
    """
    # modified from http://www.janeriksolem.net/2009/06/histogram-equalization-with-python-and.html

    # ignore zeroes
    if nozero:
        im2 = im[np.not_equal(im, 0)]
    else:
        im2 = im

    # remove extreme values
    if p:
        max = np.percentile(im2.flatten(), 100-p)
        min = np.percentile(im2.flatten(), p)
        im2[np.where(im2 > max)] = max
        im2[np.where(im2 < min)] = min

    # get image histogram
    image_histogram, bins = np.histogram(im2.flatten(), bins=nbins, density=True)
    cdf = image_histogram.cumsum()  # cumulative distribution function
    cdf = 255 * cdf / cdf[-1]  # normalize
    # use linear interpolation of cdf to find new pixel values
    image_equalized = np.interp(im.flatten(), bins[:-1], cdf)
    return image_equalized.reshape(im.shape), cdf


def read_sen2_rgb(rgbfiles, enhance=True):
    '''
    reads in 3 separate geotiff files as R G and B channels
    rgbfiles: list of three filenames including directory structure
    enhance = True: applies histogram stretching (optional)
    returns a data frame scaled to unsigned 8 bit integer values
    '''
    # make array of 8-bit unsigned integers to be memory efficient
    # open the first file with GDAL to get dimensions
    ds = gdal.Open(rgbfiles[0])
    data = ds.ReadAsArray()
    rgbdata = np.zeros([len(bands), data.shape[0], data.shape[1]], dtype=np.uint8)

    for i, thisfile in enumerate(rgbfiles):
        print('Reading data from ' + thisfile)

        # open the file with GDAL
        ds = gdal.Open(thisfile)
        data = ds.ReadAsArray()

        # only process single-band files, these have not got 3 bands
        if data.shape[0] > 3:
            # histogram stretching and keeping the values in
            #   the RGB data array as 8 bit unsigned integers
            rgbdata[i, :, :] = np.uint8(stretch(data)[0], p=2)

        ds = None
    return rgbdata


def draw_scale_bar(ax, tifproj, bars=4, length=None, location=(0.1, 0.8), linewidth=5, col='black', zorder=20):
    """
    Plot a nice scale bar with 4 subdivisions on an axis linked to the map scale.

    ax is the axes to draw the scalebar on.
    tifproj is the map projection
    bars is the number of subdivisions of the bar (black and white chunks)
    length is the length of the scalebar in km.
    location is left side of the scalebar in axis coordinates.
    (ie. 0 is the left side of the plot)
    linewidth is the thickness of the scalebar.
    color is the color of the scale bar and the text

    modified from
    https://stackoverflow.com/questions/32333870/how-can-i-show-a-km-ruler-on-a-cartopy-matplotlib-plot/35705477#35705477

    """
    # Get the limits of the axis in map coordinates
    x0, x1, y0, y1 = ax.get_extent(tifproj)

    # Set the relative position of the scale bar
    sbllx = x0 + (x1 - x0) * location[0]
    sblly = y0 + (y1 - y0) * location[1]

    # Turn the specified relative scalebar location into coordinates in metres
    sbx = x0 + (x1 - x0) * location[0]
    sby = y0 + (y1 - y0) * location[1]

    # Get the thickness of the scalebar
    thickness = (y1 - y0) / 20

    # Calculate a scale bar length if none has been given
    if not length:
        length = (x1 - x0) / 1000 / bars  # in km
        ndim = int(np.floor(np.log10(length)))  # number of digits in number
        length = round(length, -ndim)  # round to 1sf

        # Returns numbers starting with the list
        def scale_number(x):
            if str(x)[0] in ['1', '2', '5']:
                return int(x)
            else:
                return scale_number(x - 10 ** ndim)

        length = scale_number(length)

    # Generate the x coordinate for the ends of the scalebar
    bar_xs = [sbx, sbx + length * 1000 / bars]

    # Generate the y coordinate for the ends of the scalebar
    bar_ys = [sby, sby + thickness]

    # Plot the scalebar chunks
    barcol = 'white'
    for i in range(0, bars):
        # plot the chunk
        rect = patches.Rectangle((bar_xs[0], bar_ys[0]), bar_xs[1] - bar_xs[0], bar_ys[1] - bar_ys[0],
                                 linewidth=1, edgecolor='black', facecolor=barcol)
        ax.add_patch(rect)

        #        ax.plot(bar_xs, bar_ys, transform=tifproj, color=barcol, linewidth=linewidth, zorder=zorder)

        # alternate the colour
        if barcol == 'white':
            barcol = col
        else:
            barcol = 'white'
        # Generate the x,y coordinates for the number
        bar_xt = sbx + i * length * 1000 / bars
        bar_yt = sby + thickness

        # Plot the scalebar label for that chunk
        ax.text(bar_xt, bar_yt, str(round(i * length / bars)), transform=tifproj,
                horizontalalignment='center', verticalalignment='bottom',
                color=col, zorder=zorder)
        # work out the position of the next chunk of the bar
        bar_xs[0] = bar_xs[1]
        bar_xs[1] = bar_xs[1] + length * 1000 / bars
    # Generate the x coordinate for the last number
    bar_xt = sbx + length * 1000
    # Plot the last scalebar label
    t = ax.text(bar_xt, bar_yt, str(round(length)) + ' km', transform=tifproj,
                horizontalalignment='center', verticalalignment='bottom',
                color=col, zorder=zorder)


def map_it(rgbdata, imgproj, imgextent, shapefile, mapfile='map.jpg',
           maptitle='', figsizex=8, figsizey=8, zoom=1, xoffset=0, yoffset=0):
    '''
    New map_it function with scale bar located below the map but inside the enlarged map area
    This version creates different axes objects for the map, the location map and the legend.

    rgbdata = numpy array of the red, green and blue channels, made by read_sen2rgb
    imgproj = map projection of the tiff files from which the rgbdata originate
    imgextent = extent of the satellite image in map coordinates
    shapefile = shapefile name to be plotted on top of the map
    plotfile = output filename for the map plot
    plottitle = text to be written above the map
    figsizex = width of the figure in inches
    figsizey = height of the figure in inches
    zoom = zoom factor
    xoffset = offset in x direction in pixels
    yoffset = offset in x direction in pixels

    ax1 is the axes object for the main map area
    ax2 is the axes object for the location overview map in the bottom left corner
    ax3 is the axes object for the entire figure area
    ax4 is the axes object for the north arrow
    ax5 is the axes object for the map legend
    ax6 is the axes object for the map title
    mapextent = extent of the map to be plotted in map coordinates
    shpproj = map projection of the shapefile

    '''

    # work out the map extent based on the image extent plus a margin
    width = (imgextent[1] - imgextent[0]) * zoom  # work out the width and height of the zoom image
    height = (imgextent[3] - imgextent[2]) * zoom
    cx = (imgextent[0] + imgextent[1]) / 2 + xoffset  # calculate centre point positions
    cy = (imgextent[2] + imgextent[3]) / 2 + yoffset
    mapextent = (cx - width / 2, cx + width / 2, cy - height / 2, cy + height / 2)  # create a new tuple 'mapextent'

    # get shapefile projection from the file
    # get driver to read a shapefile and open it
    driver = ogr.GetDriverByName('ESRI Shapefile')
    dataSource = driver.Open(shapefile, 0)
    if dataSource is None:
        sys.exit('Could not open ' + shapefile)  # exit with an error code
    # get the layer from the shapefile
    layer = dataSource.GetLayer()

    # get the projection information and convert to wkt
    projsr = layer.GetSpatialRef()
    #print(projsr)
    projwkt = projsr.ExportToWkt()
    #print(projwkt)
    projosr = osr.SpatialReference()
    # convert wkt projection to Cartopy projection
    projosr.ImportFromWkt(projwkt)
    #print(projosr)
    projcs = projosr.GetAuthorityCode('PROJCS')
    if projcs == None:
        print("No EPSG code found in shapefile. Using EPSG 4326 instead. Make sure the .prj file contains AUTHORITY={CODE}.")
        projcs = 4326 # if no EPSG code given, set to geojson default
    print(projcs)
    if projcs == 4326:
        shapeproj = ccrs.PlateCarree()
    else:
        shapeproj = ccrs.epsg(projcs)   # Returns the projection which corresponds to the given EPSG code.
                                        # The EPSG code must correspond to a “projected coordinate system”,
                                        # so EPSG codes such as 4326 (WGS-84) which define a “geodetic
                                        # coordinate system” will not work.
    print("\nShapefile projection:")
    print(shapeproj)

    # make the figure
    fig = plt.figure(figsize=(figsizex, figsizey))

    # ---------------------- Surrounding frame ----------------------
    # set up frame full height, full width of figure, this must be called first
    left = -0.01
    bottom = -0.01
    width = 1.02
    height = 1.02
    rect = [left, bottom, width, height]
    ax3 = plt.axes(rect)

    # turn on the spines we want
    blank_axes(ax3)
    ax3.spines['right'].set_visible(False)
    ax3.spines['top'].set_visible(False)
    ax3.spines['bottom'].set_visible(False)
    ax3.spines['left'].set_visible(False)

    # add copyright statement and production date in the bottom left corner
    ax3.text(0.03, 0.03, copyright +
             'Map generated at ' + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
             fontsize=9)

    # ---------------------- Main map ----------------------
    # set up main map almost full height (allow room for title), to the right of the figure
    left = 0.3
    bottom = 0.01
    width = 0.69
    height = 0.87

    rect = [left, bottom, width, height]
    ax1 = plt.axes(rect, projection=imgproj, )

    # add 10% margin below the main map area of the image
    extent1 = (mapextent[0], mapextent[1],
               mapextent[2] - 0.1 * (mapextent[3] - mapextent[2]), mapextent[3])
    ax1.set_extent(extent1, crs=imgproj)

    RIVERS_10m = cartopy.feature.NaturalEarthFeature('physical', 'rivers_lake_centerlines', '10m',
                                                     edgecolor='blue',facecolor='none')
    BORDERS2_10m = cartopy.feature.NaturalEarthFeature('cultural', 'admin_1_states_provinces',
                                                       '10m', edgecolor='red', facecolor='none',
                                                       linestyle='-')

    ax1.add_feature(RIVERS_10m, zorder=1.2)
    ax1.add_feature(cartopy.feature.COASTLINE, edgecolor='gray', color='none', zorder=1.2)
    ax1.add_feature(BORDERS2_10m, zorder=1.2)
    ax1.stock_img()

    # work out gridline positions
    xticks, yticks = get_gridlines(mapextent[0], mapextent[1], mapextent[2], mapextent[3], nticks=6)
    gl = ax1.gridlines(crs=imgproj, xlocs=xticks, ylocs=yticks, linestyle='--', color='grey',
                       alpha=1, linewidth=1, zorder=1.3)
    ax1.set_xticks(xticks[1:-1], crs=imgproj)
    ax1.set_yticks(yticks[1:-1], crs=imgproj)

    # set axis tick mark parameters
    ax1.tick_params(bottom=False, top=True, left=True, right=False,
                    labelbottom=False, labeltop=True, labelleft=True, labelright=False)
    # N.B. note that zorder of axis ticks is reset to he default of 2.5 when the plot is drawn. This is a known bug.

    # rotate x axis labels
    ax1.tick_params(axis='x', labelrotation=90)

    # show the data from the RGB image
    img = ax1.imshow(rgbdata[:3, :, :].transpose((1, 2, 0)),
                     extent=imgextent, origin='upper', zorder=1)

    #  read shapefile and plot it onto the tiff image map
    shape_feature = ShapelyFeature(Reader(shapefile).geometries(), crs=shapeproj,
                                   edgecolor='yellow', linewidth=2,
                                   facecolor='none')
    # higher zorder means that the shapefile is plotted over the image
    ax1.add_feature(shape_feature, zorder=1.2)

    # ------------------------scale bar ----------------------------
    # adapted from https://stackoverflow.com/questions/32333870/how-can-i-show-a-km-ruler-on-a-cartopy-matplotlib-plot/35705477#35705477

    bars = 4 # plot four bar segments

    # Get the limits of the axis in map coordinates
    x0, x1, y0, y1 = ax1.get_extent(crs=imgproj) # get axes extent in map coordinates
    length = (x1 - x0) / 1000 / 3 / bars # in km    # length of scale bar segments adds up to 33% of the map width
    ndim = int(np.floor(np.log10(length)))  # number of digits in number
    length = round(length, -ndim) # round to 1sf

    # Returns numbers starting with the list
    def scale_number(x):
        if str(x)[0] in ['1', '2', '5']:
            return int(x)
        else:
            return scale_number(x - 10 ** ndim)
    length = scale_number(length)

    # relative scalebar location in map coordinates, e.g. metres
    sbx = x0 + 0.01 * (x1 - x0)
    sby = y0 + 0.04 * (y1 - y0)

    # thickness of the scalebar
    thickness = (y1 - y0) / 80

    # Generate the xy coordinates for the ends of the scalebar segment
    bar_xs = [sbx, sbx + length * 1000]
    bar_ys = [sby, sby + thickness]

    # Plot the scalebar chunks
    barcol = 'white'
    for i in range(0, bars):
        # plot the chunk
        rect = patches.Rectangle((bar_xs[0], bar_ys[0]), bar_xs[1] - bar_xs[0], bar_ys[1] - bar_ys[0],
                                 linewidth=1, edgecolor='black', facecolor=barcol, zorder=4)
        ax1.add_patch(rect)

        # alternate the colour
        if barcol == 'white':
            barcol = 'black'
        else:
            barcol = 'white'
        # Generate the x,y coordinates for the number
        bar_xt = sbx + i * length * 1000
        bar_yt = sby + thickness

        # Plot the scalebar label for that chunk
        ax1.text(bar_xt, bar_yt, str(i * length), transform=imgproj,
                 horizontalalignment='center', verticalalignment='bottom', color='black', zorder=4)

        # work out the position of the next chunk of the bar
        bar_xs[0] = bar_xs[1]
        bar_xs[1] = bar_xs[1] + length * 1000

    # Generate the x,y coordinates for the last number annotation
    bar_xt = sbx + bars * length * 1000
    bar_yt = sby + thickness

    # Plot the last scalebar label
    ax1.text(bar_xt, bar_yt, str(length * bars), transform=imgproj,
             horizontalalignment='center', verticalalignment='bottom', color='black', zorder=4)

    # work out xy coordinates for the position of the unit annotation
    bar_xt = sbx + length * bars * 500
    bar_yt = sby - thickness * 3
    # add the text annotation below the scalebar
    t = ax1.text(bar_xt, bar_yt, 'km', transform=imgproj,
                 horizontalalignment='center', verticalalignment='bottom', color='black', zorder=4)

    # do not draw the bounding box around the scale bar area. This seems to be the only way to make this work.
    #   there is a bug in Cartopy that always draws the box.
    ax1.outline_patch.set_visible(False)
    # remove the facecolor of the geoAxes
    ax1.background_patch.set_visible(False)
    # plot a white rectangle underneath the scale bar to blank out the background image over the bottom map extension
    rect = patches.Rectangle((x0, y0), x1 - x0, (y1 - y0) * 0.1, linewidth=1,
                             edgecolor='white', facecolor='white', zorder=3)
    ax1.add_patch(rect)

    # ---------------------------------Overview Location Map ------------------------
    # define where it should go, i.e. bottom left of the figure area
    left = 0.03
    bottom = 0.1
    width = 0.17
    height = 0.2
    rect = [left, bottom, width, height]

    # define the extent of the overview map in map coordinates
    #   get the map extent in latitude and longitude
    extll = ax1.get_extent(crs=ccrs.PlateCarree())
    margin = 5  # add n times the map extent
    mapw = extll[1] - extll[0] # map width
    maph = extll[3] - extll[2] # map height

    left2 = extll[0] - mapw * margin
    right2 = extll[1] + mapw * margin
    bottom2 = extll[2] - maph * margin
    top2 = extll[3] + maph * margin
    extent2 = (left2, right2, bottom2, top2)

    ax2 = plt.axes(rect, projection=ccrs.PlateCarree(), )
    ax2.set_extent(extent2, crs=ccrs.PlateCarree())
    #  ax2.set_global()  will show the whole world as context

    ax2.coastlines(resolution='110m', color='grey', zorder=3.5)
    ax2.add_feature(cfeature.LAND, color='dimgrey', zorder=1.1)
    ax2.add_feature(cfeature.BORDERS, edgecolor='red', linestyle='-', zorder=3)
    ax2.add_feature(cfeature.OCEAN, zorder=2)

    # overlay shapefile
    shape_feature = ShapelyFeature(Reader(shapefile).geometries(), crs=shapeproj,
                                   edgecolor='yellow', linewidth=1,
                                   facecolor='none')
    ax2.add_feature(shape_feature, zorder=4)

    ax2.gridlines(zorder=3)

    # add location box of the main map
    box_x = [x0, x1, x1, x0, x0]
    box_y = [y0, y0, y1, y1, y0]
    plt.plot(box_x, box_y, color='black', transform=imgproj, linewidth=1, zorder=6)

    # -------------------------------- Title -----------------------------
    # set up map title at top right of figure
    left = 0.2
    bottom = 0.95
    width = 0.8
    height = 0.04
    rect = [left, bottom, width, height]
    ax6 = plt.axes(rect)
    ax6.text(0.5, 0.0, maptitle, ha='center', fontsize=11, fontweight='bold')
    blank_axes(ax6)

    # ---------------------------------North Arrow  ----------------------------
    #
    left = 0.03
    bottom = 0.35
    width = 0.1
    height = 0.1
    rect = [left, bottom, width, height]
    ax4 = plt.axes(rect)

    # add a graphics file with a North Arrow
    compassrose = im.imread(rosepath + 'compassrose.jpg')
    img = ax4.imshow(compassrose, zorder=4) #origin='upper'

    # need a font that support enough Unicode to draw up arrow. need space after Unicode to allow wide char to be drawm?
    #ax4.text(0.5, 0.0, r'$\uparrow N$', ha='center', fontsize=30, family='sans-serif', rotation=0)
    blank_axes(ax4)

    # ------------------------------------  Legend -------------------------------------
    # legends can be quite long, so set near top of map
    left = 0.03
    bottom = 0.49
    width = 0.17
    height = 0.4
    rect = [left, bottom, width, height]
    ax5 = plt.axes(rect)
    blank_axes(ax5)

    # create an array of color patches and associated names for drawing in a legend
    # colors are the predefined colors for cartopy features (only for example, Cartopy names are unusual)
    colors = sorted(cartopy.feature.COLORS.keys())

    # handles is a list of patch handles
    handles = []
    # names is the list of corresponding labels to appear in the legend
    names = []

    # for each cartopy defined color, draw a patch, append handle to list, and append color name to names list
    for c in colors:
        patch = patches.Patch(color=cfeature.COLORS[c], label=c)
    handles.append(patch)
    names.append(c)
    # end for

    # do some example lines with colors
    river = mlines.Line2D([], [], color='blue', marker='',
                          markersize=15, label='river')
    coast = mlines.Line2D([], [], color='grey', marker='',
                          markersize=15, label='coast')
    bdy = mlines.Line2D([], [], color='red', marker='',
                        markersize=15, label='border')
    handles.append(river)
    handles.append(coast)
    handles.append(bdy)
    names.append('river')
    names.append('coast')
    names.append('border')

    # create legend
    ax5.legend(handles, names, loc='upper left')
    ax5.set_title('Legend', loc='left')

    # show the map
    fig.show()

    # save it to a file
    fig.savefig(mapfile)
    plt.close(fig)

#############################################################################
# MAIN
#############################################################################

# go to working directory
os.chdir(wd)

# make a 'maps' directory (if it does not exist yet) for output files
if not os.path.exists(mapdir):
    print("Creating directory: ", mapdir)
    os.mkdir(mapdir)

# get Sentinel L2A scene list from data directory
allscenes = [f for f in listdir(datadir) if isdir(join(datadir, f))]
allscenes = sorted(allscenes)
print('\nSentinel-2 directory: ' + datadir)
print('\nList of Sentinel-2 scenes:')
for scene in allscenes:
    if not(scene.endswith('.SAFE')):
        allscenes.remove(scene)  # remove all directory names except SAFE files
    else:
        print(scene)
print('\n')

if len(allscenes) > 0:
    for x in range(len(allscenes)):
        print("Caracas")
        scenedir = datadir + allscenes[x] + "/"
        print("Reading scene", x + 1, ":", scenedir)
        os.chdir(scenedir) # set working directory to the Sentinel scene subdirectory
        # to get the spatial footprint of the scene from the metadata file:
        # get the list of filenames ending in .xml, but exclude 'INSPIRE.xml'
        xmlfiles = [f for f in os.listdir(scenedir) if f.endswith('.xml') & (1 - f.startswith('INSPIRE'))]
        # print('Reading footprint from ' + xmlfiles[0])
        with open(xmlfiles[0], errors='ignore') as f: # use the first .xml file in the directory
            content = f.readlines()
        content = [x.strip() for x in content] # remove whitespace characters like `\n` at the end of each line
        footprint = [x for x in content if x.startswith('<EXT_POS_LIST>')] # find the footprint in the metadata
        footprint = footprint[0].split(" ") # the first element is a string, extract and split it
        footprint[0] = footprint[0].split(">")[1] #   and split off the metadata text
        footprint = footprint[:-1] #   and remove the metadata text at the end of the list
        footprint = [float(s) for s in footprint] # convert the string list to floats
        footprinty = footprint[0::2]  # list slicing to separate latitudes: list[start:stop:step]
        footprintx = footprint[1::2]  # list slicing to separate longitudes: list[start:stop:step]
        os.chdir(datadir + allscenes[x] + "/" + "GRANULE" + "/")
        sdir = listdir()[0]  # only one subdirectory expected in this directory
        imgdir = datadir + allscenes[x] + "/" + "GRANULE" + "/" + sdir + "/" + "IMG_DATA/R10m/"
        os.chdir(imgdir) # go to the image data subdirectory
        sbands = sorted([f for f in os.listdir(imgdir) if f.endswith('.jp2')]) # get the list of jpeg filenames
        print('Bands in granule directory: ')
        for band in sbands:
            print(band)
        print('Retain bands with file name pattern matching:')
        for band in bands:
            print(band)
        rgbbands = []
        for band in bands:
            goodband = [x for x in sbands if band in x]
            print(goodband)
            rgbbands.append(goodband)
        print('Band files for map making:')
        for band in rgbbands:
            print(band)
        nbands = len(rgbbands)
        if not nbands == 3:
            print("Error: Number of bands must be 3 for RGB.")
            break
        for i, iband in enumerate(rgbbands):
            print("Reading data from band " + str(i) + ": " + iband[0])
            bandx = gdal.Open(iband[0], gdal.GA_Update) # open a band
            data = bandx.ReadAsArray()
            print("Band data shape: ")
            print(data.shape)
            if i == 0:
                ncols = bandx.RasterXSize
                nrows = bandx.RasterYSize
                geotrans = bandx.GetGeoTransform()
                proj = bandx.GetProjection()
                inproj = osr.SpatialReference()
                inproj.ImportFromWkt(proj)
                ulx = geotrans[0]  # Upper Left corner coordinate in x
                uly = geotrans[3]  # Upper Left corner coordinate in y
                pixelWidth = geotrans[1]  # pixel spacing in map units in x
                pixelHeight = geotrans[5]  # (negative) pixel spacing in y
                projcs = inproj.GetAuthorityCode('PROJCS')
                projection = ccrs.epsg(projcs)
                extent = (geotrans[0], geotrans[0] + ncols * geotrans[1], geotrans[3] + nrows * geotrans[5], geotrans[3])
                rgbdata = np.zeros([nbands, data.shape[0], data.shape[1]],
                               dtype=np.uint8)  # recepticle for stretched RGB pixel values
            p=2
            print("Histogram stretching of band " + str(i) + " using p=" + str(p))
            rgbdata[i, :, :] = np.uint8(stretch(data)[0], p=p) # histogram stretching and converting to 8 bit unsigned integers
            bandx = None # close GDAL file
        #rgbdata[:,:,:] = np.uint8(stretch(rgbdata)[0]) # histogram stretching and converting to 8 bit unsigned integers

        # plot the image as RGB on a cartographic map
        # Overview map: make a map plot of the tiff file in the image projection
        mytitle = allscenes[x].split('.')[0]
        mapfile = mapdir + mytitle + '.jpg'
        print('   shapefile = ' + shapefile)
        print('   output map file = ' + mapfile)
        zoom = 1
        xoffset = 0
        yoffset = 0
        map_it(rgbdata, imgproj=projection, imgextent=extent, shapefile=shapefile,
               mapfile=mapfile, maptitle=mytitle,
               zoom=zoom, xoffset=xoffset, yoffset=yoffset)

'''
# Zoom out, i.e. zoom factor greater than 1
nfiles, mapfiles = geotif2maps(tiffroot, shapedir + shapefile, plotdir, bands=[5, 4, 3],
                               id='map2', zoom=2, xoffset=0, yoffset=0)
print('Made map files:')
for f in mapfiles: print(f)

# Zoom in to the centre, i.e. zoom factor smaller than 1
nfiles, mapfiles = geotif2maps(tiffroot, shapedir + shapefile, plotdir, bands=[5, 4, 3],
                               id='map3', zoom=1/4, xoffset=0, yoffset=0)
print('Made map files:')
for f in mapfiles: print(f)

# Zoom in to a position relative to the scene centre
nfiles, mapfiles = geotif2maps(tiffroot, shapedir + shapefile, plotdir, bands=[5, 4, 3],
                               id='map4', zoom=1/4, xoffset=round(-109800*0.25), yoffset=round(-109800*0.125))
print('Made map files:')
for f in mapfiles: print(f)

# Zoom in more
nfiles, mapfiles = geotif2maps(tiffroot, shapedir + shapefile, plotdir, bands=[5, 4, 3],
                               id='map5', zoom=1/8, xoffset=round(-109800*0.25), yoffset=round(-109800*0.1))
print('Made map files:')
for f in mapfiles: print(f)

# Zoom in even more
nfiles, mapfiles = geotif2maps(tiffroot, shapedir + shapefile, plotdir, bands=[5, 4, 3],
                               id='map6', zoom=1/16, xoffset=round(-109800*0.25), yoffset=round(-109800*0.1))
print('Made map files:')
for f in mapfiles: print(f)

# Zoom in even more
nfiles, mapfiles = geotif2maps(tiffroot, shapedir + shapefile, plotdir, bands=[5, 4, 3],
                               id='map7', zoom=1/256, xoffset=round(-109800*0.311), yoffset=round(-109800*0.134))
print('Made map files:')
for f in mapfiles: print(f)
'''
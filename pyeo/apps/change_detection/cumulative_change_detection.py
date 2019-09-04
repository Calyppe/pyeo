"""
This script takes a list of classification and calculates the mode / majority class
(i.e. most frequent class in time).
"""

import os
from glob import glob
import pyeo.core as pyeo
import numpy as np
import argparse
from scipy import stats
import matplotlib.pyplot as plt
import gdal, ogr


def write_geotiff(fname, data, geo_transform, projection, n_layers='Null'):
    """Create a GeoTIFF file with the given data. Currently only work for 1 band images.

    ToDo: Add function to automatically write out multi-band images.

    :param fname: (string) representing the output filename.
    :param data: numpy array.
    :param geo_transform:
    :param projection:
    :param int n_bnds: Number of bands in output raster. Defaults to single band images.
    """
    driver = gdal.GetDriverByName('GTiff')
    data_2d = np.delete(data,0)
    if n_layers is 'Null':
        n_bnds, rows, cols = data.shape
        dataset = driver.Create(fname, cols, rows, n_bnds, gdal.GDT_Byte)
        dataset.SetGeoTransform(geo_transform)
        dataset.SetProjection(projection)
        dataset.GetRasterBand(1).WriteArray(data[0])
        dataset = None  # Closing the file
    else:
        n_bnds, rows, cols = data.shape
        n_bnds = n_layers
        dataset = driver.Create(fname, cols, rows, n_bnds, gdal.GDT_Byte)
        dataset.SetGeoTransform(geo_transform)
        dataset.SetProjection(projection)
        dataset.GetRasterBand(1).WriteArray(data[0])
        dataset = None  # Closing the file


def crop_to_aoi(in_rst, rst_out_fn, aoi, n_bands):
    """
    These functions are written after Chris Garrard's book "Geoprocessing with Python" published by Manning Shelter
    Island, NY, in 2016.

    Thanks goes to the opensource community and to all authors who worked on the packages used in this work as well as
    to my dear colleagues who were/are key in problem solving and advancing this project.

    :param str in_rst: Path and name of raster to be processed.
    :param str rst_out_fn: Name of raster subset to be written to disk.
    :param str aoi: Path to shapefile of area of interest.
    :param int n_bands: Number of bands. Integer.
    """
    # set corner coordinates of aoi to subset image from
    shp = ogr.Open(aoi)
    shp_lyr = shp.GetLayer(0)
    extent = shp_lyr.GetExtent()

    aoi_ulx, aoi_uly = extent[0], extent[3]
    aoi_lrx, aoi_lry = extent[1], extent[2]

    in_ds = gdal.Open(in_rst)
    in_gt = in_ds.GetGeoTransform()

    # TODO: Get number of bands (nBands) from input raster to pass on to driver.Create()

    inv_gt = gdal.InvGeoTransform(in_gt)
    if gdal.VersionInfo()[0] == '1':
        if inv_gt[0] == 1:
            inv_gt = inv_gt[1]
        else:
            raise RunTimeError('Inverse geotransform failed')
    elif inv_gt is None:
        raise RunTimeError('Inverse geotransform failed')

    # Compute upper left and lower right offsets
    offsets_ul = gdal.ApplyGeoTransform(inv_gt, aoi_ulx, aoi_uly)
    offsets_lr = gdal.ApplyGeoTransform(inv_gt, aoi_lrx, aoi_lry)
    off_ulx, off_uly = map(int, offsets_ul)
    off_lrx, off_lry = map(int, offsets_lr)

    # Compute number of rows and columns to extract
    rows = off_lry - off_uly
    columns = off_lrx - off_ulx

    gtiff_driver = gdal.GetDriverByName('GTiff')
    out_ds = gtiff_driver.Create(rst_out_fn, columns, rows, n_bands)
    out_ds.SetProjection(in_ds.GetProjection())

    # Put new origin coordinates in geotransform
    subset_ulx, subset_uly = gdal.ApplyGeoTransform(in_gt, off_ulx, off_uly)
    out_gt = list(in_gt)
    out_gt[0] = subset_ulx
    out_gt[3] = subset_uly
    out_ds.SetGeoTransform(out_gt)

    if n_bands is 1:
        i = 1
        in_band = in_ds.GetRasterBand(i)
        out_band = out_ds.GetRasterBand(i)
        # Read in data using computed values
        data = in_band.ReadAsArray(off_ulx, off_uly, columns, rows)
        # Write out data starting at the origin
        out_band.WriteArray(data)

        del out_ds

    else:
        for i in range(1, n_bands):
            in_band = in_ds.GetRasterBand(i)
            out_band = out_ds.GetRasterBand(i)
            # Read in data using computed values
            data = in_band.ReadAsArray(off_ulx, off_uly, columns, rows)
            # Write out data starting at the origin
            out_band.WriteArray(data)

        del out_ds


if __name__=='__main__':
    # Reading in config file
    parser = argparse.ArgumentParser(description='Calculate the majority class (i.e. mode) using a '
                                                 'list of classifications in raster format')
    parser.add_argument("directory_of_classifications")
    # parser.add_argument("model")
    # parser.add_argument("output")
    # parser.add_argument("-l", "--log_path", default=os.path.join(os.getcwd(), ".log"))
    # parser.add_argument("-c", "--chunks", default=16)
    # parser.add_argument("-m", "--mask", action="store_true")
    args = parser.parse_args()

    # log = pyeo.init_log(args.log_path)
    directory_of_classifications = '/media/ubuntu/data_archive/F2020/Kenya/outputs/classifications/' \
                                   'kwale/cumulative_change'

    aoi_pth = '/media/ubuntu/data_archive/F2020/Kenya/vector/Kwale_extent_T37MER_new.shp'

    output_stk_pth = '/media/ubuntu/data_archive/F2020/Kenya/outputs/classifications/kwale/' \
                     'two_image_comparison/' \
                     'test_stkALLimgs.tif'

    output_stk_pth_clipped = '/media/ubuntu/data_archive/F2020/Kenya/outputs/classifications/kwale/' \
                             'two_image_comparison/' \
                             'test_stkALLimgs_clipped.tif'

    cls_mode_output = '/media/ubuntu/data_archive/F2020/Kenya/outputs/classifications/kwale/' \
                     'two_image_comparison/' \
                     'test_stkALLmode.tif'

    # get list of input raster
    rst_list = glob(os.path.join(directory_of_classifications, '*T37MER*.tif'))

    # ToDo: Clip raster to window area extent / outline.
    for id, rst in enumerate(rst_list):
        # ToDo: extract name from list and add clipped at end.
        pth_fn, fmt = rst.split(sep='.')
        out_fn = pth_fn + '_clipped' + '.' + fmt
        crop_to_aoi(in_rst=rst, rst_out_fn=out_fn, aoi=aoi_pth, n_bands=1)

    clipped_rst_list = glob(os.path.join(directory_of_classifications, '*T37MER*_clipped.tif'))

    # stack every 4th or 5th available classification
    pyeo.stack_images(clipped_rst_list, output_stk_pth)

    # ToDo: Check what frequency (i.e. when class is confirmed) gives best results.
    # -> Maybe best result two image change map.

    # read in output raster as array
    rst_stk = pyeo.raster_to_array(output_stk_pth)

    # get majority class (i.e. mode of classes) using scipy
    arr_mode, arr_count = stats.mode(rst_stk)

    # omitting NaNs
    rst_arr_2float = rst_stk.astype('float')
    # can't convert integers to nan => mode calculation takes more time
    rst_arr_2float[rst_arr_2float == 0] = np.nan
    # calc mode, omitting NaNs. Returns a masked array.
    arr_mode_nan, arr_count_nan = stats.mode(rst_arr_2float, nan_policy='omit')
    # convert nan back to 0s

    # convert back to integer 32
    arr_mode_nan_int32 = arr_mode_nan.astype('int32')

    # get mask of masked_array or return full boolean array of False
    msk = np.ma.getmask(arr_mode_nan_int32)

    # ToDo: convert fill_value=-2147483648 back to 0s.


    # get geo_transform and projection info from original reference raster
    ds = gdal.Open(rst_list[0])
    gt = ds.GetGeoTransform()
    prj = ds.GetProjection()

    # write out as geotiff
    write_geotiff(fname=cls_mode_output, data=arr_mode_nan_int32, geo_transform=gt,
                  projection=prj)





    # Kwale_T37MER.shp extent:
    # 507413.7199224829673767,9475385.7046883776783943 :
    # 577015.2497345693409443,9557457.9651870708912611
    # Kwale.shp
    # 507413.7199224829673767,9475385.7046883776783943 :
    # 577015.2497345693409443,9557457.9651870708912611



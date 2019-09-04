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
import gdal


def write_geotiff(fname, data, geo_transform, projection, n_layers='Null'):
    """Create a GeoTIFF file with the given data.

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

    output_stk_pth = '/media/ubuntu/data_archive/F2020/Kenya/outputs/classifications/kwale/' \
                     'two_image_comparison/' \
                     'test_stkALLimgs.tif'

    cls_mode_output = '/media/ubuntu/data_archive/F2020/Kenya/outputs/classifications/kwale/' \
                     'two_image_comparison/' \
                     'test_stkALLmode.tif'

    # get list of input raster
    rst_list = glob(os.path.join(directory_of_classifications, '*.tif'))

    # stack every 4th or 5th available classification
    pyeo.stack_images(rst_list, output_stk_pth)

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

    # get geo_transform and projection info from original reference raster
    ds = gdal.Open(rst_list[0])
    gt = ds.GetGeoTransform()
    prj = ds.GetProjection()

    # write out as geotiff
    write_geotiff(fname=cls_mode_output, data=arr_mode_nan_int32, geo_transform=gt,
                  projection=prj)









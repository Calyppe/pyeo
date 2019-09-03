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
                     'test_stk4imgs.tif'

    # get list of input raster
    rst_list = glob(os.path.join(directory_of_classifications, '*.tif'))

    # stack every 4th or 5th available classification
    pyeo.stack_images(rst_list[0:4], output_stk_pth)
    # read in output raster as array
    rst_stk = pyeo.raster_to_array(output_stk_pth)
    # get majority class (i.e. mode of classes) using scipy
    arr_mode, arr_count = stats.mode(rst_stk)
    plt.show(arr_mode)









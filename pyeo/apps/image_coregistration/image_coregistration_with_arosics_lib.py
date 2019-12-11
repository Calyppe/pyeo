"""
A suit of image co-registration tests using AROSICS library developed by Scheffler et al. (2017):

AROSICS: An Automated and Robust Open-Source Image Co-Registration Software for Multi-Sensor Satellite Data.

Paper link: https://www.mdpi.com/2072-4292/9/7/676

PyPi link: https://pypi.org/project/arosics/

PyDocs link: http://danschef.gitext.gfz-potsdam.de/arosics/doc/
"""
import os
from geoarray import GeoArray  # GeoArray allows to calculate shifts without disk access .
from arosics import COREG
from arosics import DESHIFTER

import pyeo.raster_manipulation as rst

# First library test
# calculate spatial shifts - with input data on disk
im_reference = '/media/ubuntu/data_archive/F2020/Colombia/raster/cartagena/images/fmas_sen2cor_combined/' \
               'S2A_MSIL2A_20180129T153111_N0206_R025_T18NWF_20180129T185058.tif'
im_target = '/media/ubuntu/data_archive/F2020/Colombia/raster/cartagena/images/' \
            'fmas_sen2cor_combined/S2A_MSIL2A_20180208T152641_N0206_R025_T18NWF_20180208T185301.tif'

target_out = '/media/ubuntu/data_archive/F2020/Colombia/raster/cartagena/images/fmas_sen2cor_combined/coreg/' \
             'S2A_MSIL2A_20180208T152641_N0206_R025_T18NWF_20180208T185301_PYEO_test.tif'

CR = COREG(im_reference, im_target, path_out=target_out, fmt_out='GTIFF', ws=(256, 256))
CR.calculate_spatial_shifts()
DESHIFTER(im_target, CR.coreg_info, path_out=target_out).correct_shifts()

# Test PYEO implementation.
cr, img_pth = rst.global_coreg(im_reference, im_target, out_fn=target_out)




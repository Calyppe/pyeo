"""A small set of functions for producing validation points from maps"""

import numpy as np
import gdal
import random
import ogr, osr
from pyeo import core
import logging
import datetime

gdal.UseExceptions()


def produce_stratifed_validation_points(map_path, n_points, out_path, no_data=None, seed=None):
    """Produces a set of stratified validation points from map_path"""
    log = logging.getLogger(__name__)
    log.info("Producing random sampling of {} with {} points.".format(map_path, n_points))
    map = gdal.Open(map_path)
    gt = map.GetGeoTransform()
    proj = map.GetProjection()
    map = None
    point_list = stratified_random_sample(map_path, n_points, no_data, seed)
    save_point_list_to_shapefile(point_list, out_path, gt, proj)
    log.info("Complete. Output saved at {}.".format(out_path))


def save_point_list_to_shapefile(point_list, out_path, geotransform, projection_wkt):
    """Saves a list of points to a shapefile at out_path. Need the gt and projection of the raster.
    GT is needed to move each point to the centre of the pixel."""
    log = logging.getLogger(__name__)
    log.info("Saving point list to shapefile")
    log.debug("GT: {}\nProjection: {}".format(geotransform, projection_wkt))
    driver = ogr.GetDriverByName("ESRI Shapefile")
    data_source = driver.CreateDataSource(out_path)
    srs = osr.SpatialReference()
    srs.ImportFromWkt(projection_wkt)
    layer = data_source.CreateLayer("validation_points", srs, ogr.wkbPoint)
    for point in point_list:
        feature = ogr.Feature(layer.GetLayerDefn())
        coord = core.pixel_to_point_coordinates(point, geotransform)
        offset = geotransform[1]/2   # Adds half a pixel offset so points end up in the center of pixels
        wkt = "POINT({} {})".format(coord[0]+offset, coord[1]+offset)
        new_point = ogr.CreateGeometryFromWkt(wkt)
        feature.SetGeometry(new_point)
        layer.CreateFeature(feature)
        feature = None
    layer = None
    data_source = None


def stratified_random_sample(map_path, n_points, minimum_class_pixels=None, no_data=None, seed = None):
    """Produces a stratified list of pixel coordinates. WARNING: high mem!"""
    log = logging.getLogger(__name__)
    if not seed:
        seed = datetime.datetime.now().timestamp()
    map = gdal.Open(map_path)
    map_array = map.GetVirtualMemArray()
    class_dict = build_class_dict(map_array, no_data)
    map_array = None
    n_pixels = sum(len(coord_list) for coord_list in class_dict.values())
    out_coord_list = []
    for pixel_class, coord_list in class_dict.items():
        proportion = len(coord_list)/n_pixels
        n_sample_pixels = int(np.round(proportion*n_points))
        if minimum_class_pixels:
            if n_sample_pixels < minimum_class_pixels[pixel_class]:
                n_sample_pixels = minimum_class_pixels[pixel_class]
        out_coord_list.extend(random.sample(coord_list, n_sample_pixels))
    return out_coord_list


def build_class_dict(class_array, no_data=None):
    """Returns a dict of coordinates of the following shape:
    [class, coord].
    WARNING: This will take up a LOT of memory!"""
    out_dict = {}
    it = np.nditer(class_array, flags=['multi_index'])
    while not it.finished:
        this_class = int(it.value)
        if this_class == no_data:
            it.iternext()
            continue
        if this_class in out_dict.keys():
            out_dict[this_class].append(it.multi_index)
        else:
            out_dict.update({this_class: [it.multi_index]})
        it.iternext()
    return out_dict


def create_validation_set(map_path, ):
    pass


def cal_si(ui):
    si = np.sqrt(ui * (1 - ui))
    return si


def cal_wi(n,total_n):
    wi = float(n/total_n)
    return wi


def cal_w_all(dict_pixel_numbers):
    w_dict = {}
    total_pixel = (sum(dict_pixel_numbers.values()))
    for key in dict_pixel_numbers:
        w_dict[key] = cal_wi(n=dict_pixel_numbers[key], total_n= total_pixel)
    return w_dict


def cal_n_by_prop(weight, sample_size):
    n = round(weight * sample_size)
    return n


def val_to_sd(val):
    sd = val**0.5
    return sd


def cal_val_for_overall_accruacy(weight_dict,u_dict,sample_size_dict):
    sum_val = 0
    for key in u_dict:
        val_i = (weight_dict[key] **2) * u_dict[key] * (1-u_dict[key])/(sample_size_dict[key]-1)
        sum_val += val_i
    return sum_val


def cal_val_for_user_accuracy(u_i,sample_size_i):
    val_user = (u_i*(1-u_i))/(sample_size_i-1)
    return val_user


def cal_sd_for_overall_accruacy(weight_dict,u_dict,sample_size_dict):
    val_overall = cal_val_for_overall_accruacy(weight_dict=weight_dict,u_dict=u_dict,sample_size_dict=sample_size_dict)
    sd_overall = val_to_sd(val_overall)
    return sd_overall


def cal_sd_for_user_accuracy(u_i,sample_size_i):
    val_user = cal_val_for_user_accuracy(u_i=u_i,sample_size_i=sample_size_i)
    sd_user = val_to_sd(val_user)
    return sd_user


def cal_total_sample_size(se_expected_overall, user_accuracy, pixel_numbers, type ='simple'):
    """
    Calculates the total sample size between all classes
    Parameters
    ----------
    se_expected_overall: The standard error
    user_accuracy
    pixel_numbers
    type

    Returns
    -------
    The minimum sample size
    """
    total_pixel = (sum(pixel_numbers.values()))
    if type == 'simple':
        weighted_U_sum = 0
        # weight are equal between different classes
        for key in user_accuracy:
            S_i = cal_si(user_accuracy[key])
            Wi = cal_wi(n= pixel_numbers[key], total_n=total_pixel)# proportion of each class
            weighted_U_sum += S_i*Wi
        n = (weighted_U_sum/se_expected_overall)**2
    elif type == 'full':
        weighted_U_sum2 = 0
        weighted_U_sum = 0
        # weight are equal between different classes
        for key in user_accuracy:
            S_i = cal_si(user_accuracy[key])
            Wi = cal_wi(n= pixel_numbers[key], total_n=total_pixel)  # proportion of each class
            weighted_U_sum2 += S_i * Wi
            weighted_U_sum += (S_i ** 2) * Wi
        up = (weighted_U_sum2) ** 2
        bottom_right = (1 / total_pixel) * weighted_U_sum

        n = (up / (se_expected_overall ** 2 + bottom_right))
    print('suggested total sample size are:' + str(n))
    return n


def calc_minimum_n(expected_accuracy, variance_tolerance):
    """
    Calculates the rminimum number of points required to achieve the specified accuracy
    Parameters
    ----------
    expected_accuracy: Between 0 and 1
    variance_tolerance:

    Returns
    -------

    """
    n = expected_accuracy * (1-expected_accuracy) / variance_tolerance
    return n


def allocate_category_sample_sizes(total_sample_size, user_accuracy, pixel_numbers, variance_tolerance,
                                   allocate_type='olofsson'):
    """
    Allocates a number of pixels to sample per class that will fulfil the parameters given

    Parameters
    ----------
    total_sample_size: The total number of validation points requested (from cal_total_sample_size)
    user_accuracy: Dictionary of estimated user accuracies for classes in map (between 0 and 1)
    pixel_numbers: Dictionary of total pixels for each class in user_accuracy
    variance_tolerance: Acceptable vairance between the sample accuary and the data accuracy with a certain sample size
    allocate_type: The allocation strategy to be used. Can be 'equal', 'prop' or 'olofsson'.

    Returns
    -------
    A dictionary of classes and no. pixels per class.

    """
    log = logging.getLogger(__name__)
    minimum_n = {}
    allocated_n = {}
    weight = cal_w_all(pixel_numbers)
    log.info('the weight for each class is: ')
    log.info(weight)
    log.info('-----------------')
    log.info('the minimum sampling number for : ')
    for key in user_accuracy:
        minimum_n_i = calc_minimum_n(expected_accuracy=user_accuracy[key], variance_tolerance=variance_tolerance)
        log.info('      ' + key + ' is: ' + str(round(minimum_n_i)))
        minimum_n[key] = minimum_n_i

    if allocate_type == 'equal' or allocate_type == 'prop':
        for key in user_accuracy:
            if allocate_type == 'equal':
                n = total_sample_size/float(len(user_accuracy.keys()))
            elif allocate_type == 'prop':
                n = cal_n_by_prop(weight=weight[key], sample_size=total_sample_size)
            else:
                continue
            allocated_n[key] = n
            log.info('allocated sampling number for ' + key + ' is: ' + str(allocated_n[key]))

    elif allocate_type == 'olofsson':
        allocated_n = part_fixed_value_sampling(allocated_n,  total_sample_size, weight)
        log.info('allocated sample number under different scenario is: ')
        log.info(allocated_n)
    else:
        raise core.ForestSentinelException("Invalid allocation type: valid values are 'equal', 'prop' or 'olofsson")
    return allocated_n


def part_fixed_value_sampling(pinned_sample_numbers, class_total_sizes, total_sample_size):
    """
    Given a dictionary of classes in
    Parameters
    ----------
    class_sample_numbers: Dictionary of class label and pinned sample numbers or 'None'.
    total_sample_size
    weight

    Returns
    -------

    """
    pinned_sample_total = sum(sample_size for sample_size in pinned_sample_numbers.values() if sample_size is not None)
    total_map_size = sum(class_total_sizes.values())
    remaining_samples = total_sample_size - pinned_sample_total
    out_values = {}
    for map_class, sample_points in pinned_sample_numbers.items():
        if sample_points is not None:
            out_values.update({map_class: sample_points})
        else:
            weight = class_total_sizes[map_class]/total_map_size
            sample_points = weight*remaining_samples
            out_values.update({map_class: sample_points})
    return out_values
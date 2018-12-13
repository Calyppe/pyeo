import sys, os
sys.path.insert(0, os.path.abspath(os.path.join(__file__, '..', '..', '..', '..')))
import pyeo.core as pyeo


def model_comparison(in_tifs, in_shps, class_attribue, model_list, out_dir):
    """Takes in a set of training data (in the standard Pyeo form of shapefile-tif) and a list
    of models. Trains each model using the same training data, and produces a classificatin map
    for each."""
    for model in model_list:
        classes = []
        signatures = []
        for tif, shp in zip(in_tifs, in_shps):
            tmp_classes, tmp_sigs = pyeo.get_training_data(tif)
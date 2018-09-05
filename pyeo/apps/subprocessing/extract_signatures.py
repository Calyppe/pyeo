import submodules as pc
import csv
import os
import glob
import numpy as np

path_to_region = r"/home/ubuntu/caqueta/training_data"

image_glob = os.path.join(path_to_region, r"*.tif")
image_list = glob.glob(image_glob)
#image_list = [r"/home/ubuntu/caqueta/training_data/colombia_cartagena__20170302_20161222.tif"]


for training_image_file_path in image_list:
    training_image_folder, training_image_name = os.path.split(training_image_file_path)
    training_image_name = training_image_name[:-4]  # Strip the file extension
    shape_path = os.path.join(training_image_folder, training_image_name, training_image_name + '.shp')
    this_training_data, this_classes = pc.get_training_data(training_image_file_path, shape_path)

    sigs=np.vstack((this_classes, this_training_data.T))

    sig_out_path = r"/home/ubuntu/caqueta/training_data/signatures/{}_signatures.csv".format(training_image_name)

    with open(sig_out_path, 'w', newline='') as outfile:
        writer = csv.writer(outfile)
        writer.writerows(sigs.T)
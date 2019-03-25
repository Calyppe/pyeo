import os, sys
sys.path.insert(0, os.path.abspath(os.path.join(__file__, '..', '..', '..', '..')))
from pyeo import core as pyeo
import argparse
import configparser
import os

# Load config

parser = argparse.ArgumentParser(description='Cla')
parser.add_argument('--conf', dest='config_path', action='store', default=r'change_detection.ini',
                    help="Path to the .ini file specifying the job.")
args = parser.parse_args()

conf = configparser.ConfigParser()
conf.read(args.config_path)
api_key_path = conf['planet']['api_key']
project_root = conf['planet']['root_dir']
aoi_path = conf['planet']['aoi_path']
start_time = conf['planet']['start_date']
end_time = conf['planet']['end_date']
log_path = conf['planet']['log_path']
model_path = conf['planet']['model']



pyeo.create_file_structure(project_root)
log = pyeo.init_log(log_path)

planet_image_path = os.path.join(project_root, r"images/planet")
stacked_image_path = os.path.join(project_root, r"images/stacked")
catagorised_image_path = os.path.join(project_root, r"output/categories")
probability_image_path = os.path.join(project_root, r"output/probabilities")

api_key = api_key_path  # For now

# Download two images

pyeo.planet_query(aoi_path, start_time, end_time, out_path=planet_image_path,
                  api_key=api_key_path, item_type = "PSScene4Band" , search_name="auto", asset_type="analytic", threads=5)
#image_1_path = pyeo.download_planet_image_on_day(aoi_path, date_1, planet_image_path, api_key)
#image_2_path = pyeo.download_planet_image_on_day(aoi_path, date_2, planet_image_path, api_key)

# Save RGB


# Do classification
images_to_stack = [os.path.join(dirpath,f)
                       for dirpath, dirnames, files in os.walk(planet_image_path)
                       for f in files if f.endswith('.tif')]

stacked_image_path = pyeo.stack_images(images_to_stack, stacked_image_path)
# pyeo.classify_image(stacked_image_path, model_path=model_path)

# Do the imshow transitive closure trick

# Polygonise

# Burn

# Resample


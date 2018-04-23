from sentinelsat.sentinel import SentinelAPI
import pandas as pd


def sen2_download(products_file, conf):
    """This function downloads Sentinel-2 imagery

    :param str products_file: The path to the CSV file containing the Sentinel-2 scene IDs.
    :param str conf: String object pointing to the config.ini file which includes all required parameters to access
    the SciHub and download the data.
    """
    # TODO: read in csv file with product IDs
    products = pd.read_csv(products_file)

    # TODO: Specify download location, file resuming.
    api = SentinelAPI(conf["sen2"]["user"], conf["sen2"]["pass"], 'https://scihub.copernicus.eu/dhus')
    # TODO: assign correct column to use for product IDs. Currently assuming the the first column is correct.
    api.download_all(products, conf["data"]["out_folder"])
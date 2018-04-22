from sentinelsat.sentinel import SentinelAPI


def sen2_download(products, conf):
    """This function downloads Sentinel-2 imagery

    :param str products: The path to the CSV  containing the Sentinel-2 scene IDs.
    :param str conf: String object pointing to the config.ini file which includes all required parameters to access
    the SciHub and download the data.
    """
    # TODO: Specify download location, file resuming.
    api = SentinelAPI(conf["sen2"]["user"], conf["sen2"]["pass"], 'https://scihub.copernicus.eu/dhus')
    api.download_all(products, conf["data"]["out_folder"])
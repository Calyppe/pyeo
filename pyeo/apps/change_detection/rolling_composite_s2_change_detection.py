"""
rolling_composite_s2_change_detection
-------------------------------------
An app for providing continuous change detection. Runs the following algorithm
 Step -1: Create initial composite with last date (stored in filename?)

 Step 0: Load composite

 Step 1: Download images from last date in composite until present last date

 Step 2: Preprocess each image

 Step 3: Generate cloud mask for each image

For each preprocessed image:

    Step 4: Build stack with composite

    Step 5: Classify stack

    Step 6: Update composite with last cloud-free pixel based on cloud mask

    Step 7: Update last_date of composite

 """
import sys

import pyeo.classification
import pyeo.queries_and_downloads
import pyeo.raster_manipulation
import pyeo.filesystem_utilities


import configparser
import argparse
import os
import datetime as dt


if __name__ == "__main__":
    do_all = True

    # Reading in config file
    parser = argparse.ArgumentParser(description='Downloads, preprocesses and classifies sentinel 2 images. A directory'
                                                 'structure to contain preprocessed and downloaded files will be'
                                                 ' created at the aoi_root location specified in the config file.'
                                                 'If any of the step flags (d,p,m,a,s,c,u,r) are present, only those '
                                                 'steps will be performed - otherwise all processing steps will be '
                                                 'performed.')
    parser.add_argument(dest='config_path', action='store', default=r'change_detection.ini',
                        help="A path to a .ini file containing the specification for the job. See "
                             "pyeo/apps/change_detection/change_detection.ini for an example.")
    parser.add_argument('--start_date', dest='start_date', help="Overrides the start date in the config file. Set to "
                                                                "LATEST to get the date of the last merged accquistion")
    parser.add_argument('--end_date', dest='end_date', help="Overrides the end date in the config file. Set to TODAY"
                                                            "to get today's date")
    parser.add_argument('-b', '--build_composite', dest='build_composite', action='store_true', default=False,
                        help="If present, creates a cloud-free (ish) composite between the two dates specified in the "
                             "config file.")
    parser.add_argument("--chunks", dest="num_chunks", type=int, default=10, help="Sets the number of chunks to split "
                                                                                  "images to in ml processing")
    parser.add_argument('--download_source', default="scihub", help="Sets the download source, can be scihub"
                                                                    "(default) or aws")
    parser.add_argument('--flip_stacks', action='store_true', default=False,
                        help="If present, stasks the classification stack as new(bgr), old(bgr). Default is"
                             "old(bgr), new(bgr). For compatability with old models.")
    parser.add_argument('--download_l2_data', action='store_true', default=False,
                        help="If present, skips sen2cor and instead downloads every image in the query with"
                             "both a L1 and L2 product")
    parser.add_argument('--build_prob_image', action='store_true', default=False,
                        help="If present, build a confidence map of pixels. These tend to be large.")

    parser.add_argument('-d', '--download', dest='do_download', action='store_true', default=False,
                        help='If present, perform the query and download level 1 images.')
    parser.add_argument('-p', '--preprocess', dest='do_preprocess', action='store_true',  default=False,
                        help='If present, apply sen2cor to all .SAFE in images/L1. Stores the result in images/L2')
    parser.add_argument('-m', '--merge', dest='do_merge', action='store_true', default=False,
                        help='If present, merges the blue, green, red and IR 10m rasters in each L2 safefile'
                             ' into a single 4-band raster. This will also mask and reproject'
                             ' the image to the requested projection. Stores the result in images/merged.')
    parser.add_argument('-a', '--mask', dest='do_mask', action='store_true', default=False,
                        help="Masks are created by default - this presently does nothing.")
    parser.add_argument('-s', '--stack', dest='do_stack', action='store_true', default=False,
                        help="For each image in images/merged, stacks with the composite in composite/ that is most "
                             "recent prior to the image. Stores an 8-band geotiff in images/stacked, where bands 1-4 "
                             "are the B,G,R and IR bands of the composite and band 5-8 are the B,G,R and IR bands of"
                             "the merged image.")
    parser.add_argument('-c', '--classify', dest='do_classify', action='store_true', default=False,
                        help="For each image in images/stacked, applies the classifier given in the .ini file. Saves"
                             "the outputs in output/categories.")
    parser.add_argument('-u', '--update', dest='do_update', action='store_true', default=False,
                        help="Builds a new cloud-free composite in composite/ from the latest image and mask"
                             " in images/merged")
    parser.add_argument('-r', '--remove', dest='do_delete', action='store_true', default=False,
                        help="Not implemented. If present, removes all images in images/ to save space.")

    parser.add_argument('--skip_prob_image', dest="skip_prob_image", action="store_true", default=False,
                        help="")

    args = parser.parse_args()

    # If any processing step args are present, do not assume that we want to do all steps
    if (args.do_download or args.do_preprocess or args.do_merge or args.do_stack or args.do_classify) == True:
        do_all = False

    conf = configparser.ConfigParser(allow_no_value=True)
    conf.read(args.config_path)

    sen_user = conf['sent_2']['user']
    sen_pass = conf['sent_2']['pass']
    project_root = conf['forest_sentinel']['root_dir']
    aoi_path = conf['forest_sentinel']['aoi_path']
    start_date = conf['forest_sentinel']['start_date']
    end_date = conf['forest_sentinel']['end_date']
    log_path = conf['forest_sentinel']['log_path']
    cloud_cover = conf['forest_sentinel']['cloud_cover']
    cloud_certainty_threshold = int(conf['forest_sentinel']['cloud_certainty_threshold'])
    model_path = conf['forest_sentinel']['model']
    sen2cor_path = conf['sen2cor']['path']
    composite_start_date = conf['forest_sentinel']['composite_start']
    composite_end_date = conf['forest_sentinel']['composite_end']
    epsg = int(conf['forest_sentinel']['epsg'])

    pyeo.filesystem_utilities.create_file_structure(project_root)
    pyeo.filesystem_utilities.create_file_structure(project_root)
    log = pyeo.filesystem_utilities.init_log(log_path)
    
    try:
        l1_image_dir = os.path.join(project_root, r"images/L1")
        l2_image_dir = os.path.join(project_root, r"images/L2")
        planet_image_dir = os.path.join(project_root, r"images/planet")
        merged_image_dir = os.path.join(project_root, r"images/merged")
        stacked_image_dir = os.path.join(project_root, r"images/stacked")
        catagorised_image_dir = os.path.join(project_root, r"output/categories")
        probability_image_dir = os.path.join(project_root, r"output/probabilities")
        composite_dir = os.path.join(project_root, r"composite")
        composite_l1_image_dir = os.path.join(project_root, r"composite/L1")
        composite_l2_image_dir = os.path.join(project_root, r"composite/L2")
        composite_merged_dir = os.path.join(project_root, r"composite/merged")

        if args.skip_prob_image:
            probability_image_dir = None

        if args.start_date == "LATEST":
            # This isn't nice, but returns the yyyymmdd string of the latest classified image
            start_date = pyeo.filesystem_utilities.get_image_acquisition_time(pyeo.filesystem_utilities.sort_by_timestamp(
                [image_name for image_name in os.listdir(catagorised_image_dir) if image_name.endswith(".tif")],
                recent_first=True
            )[0]).strftime("%Y%m%d")
        elif args.start_date:
            start_date = args.start_date

        if args.end_date == "TODAY":
            end_date = dt.date.today().strftime("%Y%m%d")
        elif args.end_date:
            end_date = args.end_date

        # Download and build the initial composite. Does not do by default
        if args.build_composite:
            if args.do_download or do_all:
                log.info("Downloading for initial composite between {} and {} with cloud cover <= ()".format(
                    composite_start_date, composite_end_date, cloud_cover))
                composite_products = pyeo.queries_and_downloads.check_for_s2_data_by_date(aoi_path, composite_start_date, composite_end_date,
                                                                                          conf, cloud_cover=cloud_cover)
                if args.download_l2_data:
                    log.info("Filtering query results for matching L1 and L2 products")
                    composite_products = pyeo.queries_and_downloads.filter_non_matching_s2_data(composite_products)
                    log.info("{} products remain".format(len(composite_products)))
                else:
                    log.info("Filtering query results to L1 only")
                    composite_products = pyeo.queries_and_downloads.filter_to_l1_data(composite_products)
                pyeo.queries_and_downloads.download_s2_data(composite_products, composite_l1_image_dir, composite_l2_image_dir,
                                                            source=args.download_source, user=sen_user, passwd=sen_pass, try_scihub_on_fail=True)
            if args.do_preprocess or do_all and not args.download_l2_data:
                log.info("Preprocessing composite products")
                pyeo.raster_manipulation.atmospheric_correction(composite_l1_image_dir, composite_l2_image_dir, sen2cor_path,
                                                                delete_unprocessed_image=False)
            if args.do_merge or do_all:
                log.info("Aggregating composite layers")
                pyeo.raster_manipulation.preprocess_sen2_images(composite_l2_image_dir, composite_merged_dir, composite_l1_image_dir,
                                                                cloud_certainty_threshold, epsg=epsg, buffer_size=5)
            log.info("Building initial cloud-free composite")
            pyeo.raster_manipulation.composite_directory(composite_merged_dir, composite_dir, generate_date_images=True)

        # Query and download all images since last composite
        if args.do_download or do_all:
            products = pyeo.queries_and_downloads.check_for_s2_data_by_date(aoi_path, start_date, end_date, conf, cloud_cover=cloud_cover)
            if args.download_l2_data:
                log.info("Filtering query results for matching L1 and L2 products")
                products = pyeo.queries_and_downloads.filter_non_matching_s2_data(products)
                log.info("{} products remain".format(len(products)))
            else:
                log.info("Filtering query results to L1 only")
                products = pyeo.queries_and_downloads.filter_to_l1_data(products)
            log.info("Downloading")
            pyeo.queries_and_downloads.download_s2_data(products, l1_image_dir, l2_image_dir, args.download_source, user=sen_user, passwd=sen_pass, try_scihub_on_fail=True)

        # Atmospheric correction
        if args.do_preprocess or do_all and not args.download_l2_data:
            log.info("Applying sen2cor")
            pyeo.raster_manipulation.atmospheric_correction(l1_image_dir, l2_image_dir, sen2cor_path, delete_unprocessed_image=False)

        # Aggregating layers into single image
        if args.do_merge or do_all:
            log.info("Aggregating layers")
            pyeo.raster_manipulation.preprocess_sen2_images(l2_image_dir, merged_image_dir, l1_image_dir, cloud_certainty_threshold, epsg=epsg,
                                                            buffer_size=5)

        log.info("Finding most recent composite")
        try:
            latest_composite_name = \
                pyeo.filesystem_utilities.sort_by_timestamp(

                    [image_name for image_name in os.listdir(composite_dir) if image_name.endswith(".tif")],
                    recent_first=True
                )[0]
            latest_composite_path = os.path.join(composite_dir, latest_composite_name)
            log.info("Most recent composite at {}".format(latest_composite_path))
        except IndexError:
            log.critical("Latest composite not found. The first time you run this script, you need to include the "
                         "--build-composite flag to create a base composite to work off. If you have already done this,"
                         "check that the earliest dated image in your images/merged folder is later than the earliest"
                         " dated image in your composite/ folder.")
            sys.exit(1)

        log.info("Sorting image list")
        images = \
            pyeo.filesystem_utilities.sort_by_timestamp(
                [image_name for image_name in os.listdir(merged_image_dir) if image_name.endswith(".tif")],
                recent_first=False
            )
        if not images:
            raise FileNotFoundError("No images found in {}. Did your preprocessing complete?".format(merged_image_dir))
        log.info("Images to process: {}".format(images))

        for image in images:
            log.info("Detecting change for {}".format(image))
            new_image_path = os.path.join(merged_image_dir, image)

            # Stack with preceding composite
            if args.do_stack or do_all:
                try:
                    latest_composite_path = pyeo.filesystem_utilities.get_preceding_image_path(new_image_path, composite_dir)
                except FileNotFoundError:
                    log.warning("No preceding composite found for {}, skipping.".format(new_image_path))
                    continue
                log.info("Stacking {} with composite {}".format(new_image_path, latest_composite_path))
                new_stack_path = pyeo.raster_manipulation.stack_image_with_composite(new_image_path, latest_composite_path, stacked_image_dir,
                                                                                     invert_stack=args.flip_stacks)
            #else new_stack_path =

            # Classify with composite
            if args.do_classify or do_all:
                log.info("Classifying with composite")
                new_class_image = os.path.join(catagorised_image_dir, "class_{}".format(os.path.basename(new_stack_path)))
                if args.build_prob_image:
                    new_prob_image = os.path.join(probability_image_dir, "prob_{}".format(os.path.basename(new_stack_path)))
                else:
                    new_prob_image = None
                pyeo.classification.classify_image(new_stack_path, model_path, new_class_image, new_prob_image, num_chunks=args.num_chunks,
                                                   skip_existing=True, apply_mask=True)

            # Build new composite
            if args.do_update or do_all:
                log.info("Updating composite")
                new_composite_path = os.path.join(
                    composite_dir, "composite_{}.tif".format(
                        pyeo.filesystem_utilities.get_sen_2_image_timestamp(os.path.basename(image))))
                pyeo.raster_manipulation.composite_images_with_mask(
                    (latest_composite_path, new_image_path), new_composite_path, generate_date_image=True)
                latest_composite_path = new_composite_path

        log.info("***PROCESSING END***")
    except Exception:
        log.exception("Fatal error in rolling_s2_composite chain")

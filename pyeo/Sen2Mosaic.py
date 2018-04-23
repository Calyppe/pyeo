import os
import gdal
import math
import glob2

def mosaic(path, out_name, wildcard=None, in_format=None):
   """Creates a mosaic of multiple rasters. Input raster format is GeoTiff
   at the moment. Will be changed to more formats in the near future.

   :param str path: Path to directory of input rasters.
   :param str wildcard: A string representing the term to search input raster
   of interest. Default wildcard searches through Sentinel-2 folder structure.
   :param in_format: Format of input dataset. Default is GeoTiff.
   :param out_name: Output filename as string.
   """
   if wildcard is None:
       wildcard = '*.tif'
   else:
       wildcard = wildcard

   # calculate output extent from all inputs
   in_files = glob2.glob(os.path.join(path, wildcard))
   min_x, max_y, max_x, min_y = get_extent(in_files[0])
   for fn in in_files[1:]:
       minx, maxy, maxx, miny = get_extent(fn)
       min_x = min(min_x, minx)
       max_y = max(max_y, maxy)
       max_x = max(max_x, maxx)
       min_y = min(min_y, miny)

   # calculate dimensions
   in_ds = gdal.Open(in_files[0])
   gt = in_ds.GetGeoTransform()
   # ceil() rounds up to cloasest integer to not cut off pixels accidentally
   rows = math.ceil((max_y - min_y) / -gt[5])
   columns = math.ceil((max_x - min_x) / gt[1])

   # create output
   if in_format is None:
       in_format = 'gtiff'
   else:
       in_format = in_format

   driver = gdal.GetDriverByName(in_format)
   out_ds = driver.Create(out_name, columns, rows)
   out_ds.SetProjection(in_ds.GetProjection())
   out_band = out_ds.GetRasterBand(1)

   # calculate new geotransform
   gt = list(in_ds.GetGeoTransform())
   gt[0], gt[3] = min_x, max_y
   out_ds.SetGeotransform(gt)

   for fn in in_files:
       in_ds = gdal.Open(fn)
       # get output offsets
       trans = gdal.Transformer(in_ds, out_ds, [])
       success, xyz = trans.TransformPoint(False, 0, 0)
       x, y, z = map(int, xyz)
       # copy data
       data = in_ds.GetRasterBand(1).ReadAsArray()
       out_band.WriteArray(data, x, y)

   # delete object
   del in_ds, out_band, out_ds


def get_extent(fn):
   """Returns min_x, max_y, max_x, min_y from a raster data.

   :param fn: Raster data set or path to raster input data.
   """
   ds = gdal.open(fn)
   gt = ds.GetGeoTransform()
   return (gt[0], gt[3], gt[0] + gt [1] * ds.RasterXSize,
           gt[3] + gt[5] * ds.RasterYSize)

if "__name__" == "__main__":
    # TODO: Test mosaic function with unit-testing function.

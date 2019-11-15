"""
A suit of image coregistration tests using scikit-image's RANSAC approach.

RANSAC -> random sample consensus.

Example taken from:

https://scikit-image.org/docs/dev/auto_examples/transform/plot_matching.html

"""
import pyeo.raster_manipulation
import pyeo.filesystem_utilities

import numpy as np
from matplotlib import pyplot as plt

from skimage import data
from skimage.util import img_as_float
from skimage.feature import (corner_harris, corner_subpix, corner_peaks, plot_matches)
from skimage.transform import warp, AffineTransform
from skimage.exposure import rescale_intensity
from skimage.color import rgb2gray
from skimage.measure import ransac

import pyeo.raster_manipulation


# generate synthetic checkerboard image and add gradient for the later matching
# ToDo: Could be gdal image load?
checkerboard = img_as_float(data.checkerboard())
img_orig = np.zeros(list(checkerboard.shape) + [3])
img_orig[..., 0] = checkerboard
gradient_r, gradient_c = (np.mgrid[0:img_orig.shape[0],
                          0:img_orig.shape[1]] / float(img_orig.shape[0]))

img_orig[...,1] = gradient_r
img_orig[...,2] = gradient_c
# ToDo: Turn multi-layer image into grayscale single band image.
img_orig = rescale_intensity(img_orig)
img_orig_gray = rgb2gray(img_orig)

# warp synthetic image i.e rotate original image to have example image later on for co-registration
# Note: Not necessary for coregister_image() function.
tform = AffineTransform(scale=(0.9, 0.9), rotation=0.2, translation=(20, -10))
img_warped = warp(img_orig, tform.inverse, output_shape=(200, 200))
img_warped_gray = rgb2gray(img_warped)

# extract corners using Harris' corner measure
# ToDo: Get corner coords of reference image and image to be co-registered.
coords_orig = corner_peaks(corner_harris(img_orig_gray), threshold_rel=0.001,
                           min_distance=5)
coords_warped = corner_peaks(corner_harris(img_warped_gray), threshold_rel=0.001,
                             min_distance=5)

# determine sub-pixel corner position - i.e. more accurate positioning.
# ToDo for coregister_image() function
coords_orig_subpix = corner_subpix(img_orig_gray, coords_orig, window_size=9)
coords_warped_subpix = corner_subpix(img_warped_gray, coords_warped, window_size=9)

def gaussian_weights(window_ext, sigma=1):
    """Computes the weights for RANSAC's image registration implemented in the coregister_image() function

    :returns: g - numpy array containing the weights.
    """
    y, x = np.mgrid[-window_ext:window_ext+1, -window_ext:window_ext+1]
    g = np.zeros(y.shape, dtype=np.double)
    g[:] = np.exp(-0.5 * (x**2 / sigma**2 + y**2 / sigma**2))
    g /= 2 *np.pi * sigma * sigma
    return g


def match_corner(coord, window_ext=5):
    r, c = np.round(coord).astype(np.intp)
    window_orig = img_orig[r-window_ext:r+window_ext+1,
                  c-window_ext:c+window_ext+1, :]

    # weight pixels depending on distance to center pixel
    weights = gaussian_weights(window_ext, sigma=3)
    weights = np.dstack((weights, weights, weights))

    # compute sum of squared differences to all corners in warped image
    SSDs = []
    for cr, cc in coords_warped:
        window_warped = img_warped[cr-window_ext:cr+window_ext+1,
                        cc-window_ext:cc+window_ext+1, :]
        SSD = np.sum(weights * (window_orig - window_warped)**2)
        SSDs.append(SSD)

    # use corner with minimum SSD as correspondence
    min_idx = np.argmin(SSDs)
    return coords_warped_subpix[min_idx]


# find correspondences using simple weighted sum of squared differences
src = []
dst = []
for coord in coords_orig_subpix:
    src.append(coord)
    dst.append(match_corner(coord))
src = np.array(src)
dst = np.array(dst)

# estimate affine transformation model using all coordinates
model = AffineTransform()
model.estimate(src, dst)

# robustly estimate affine transform model with RANSAC
model_robust, inliers = ransac((src, dst), AffineTransform, min_samples=3,
                               residual_threshold=2, max_trials=100)
outliers = inliers == False

# compare "true" and estimated transform parameters
print("Ground truth:")
print(f"Scale: ({tform.scale[1]:.4f}, {tform.scale[0]:.4f}), "
      f"Translation: ({tform.translation[1]:.4f}, "
      f"{tform.translation[0]:.4f}), "
      f"Rotation: ({-tform.rotation:.4f}")
print("Affine transform:")
print(f"Scale: ({model.scale[1]:.4f}, {model.scale[0]:.4f}), "
      f"Translation: ({model.translation[1]:.4f}, "
      f"{model.translation[0]:.4f}), "
      f"Rotation: ({-model.rotation:.4f}")
print("RANSAC:")
print(f"Scale: ({model_robust.scale[1]:.4f}, {model_robust.scale[0]:.4f}), "
      f"Translation: ({model_robust.translation[1]:.4f}, "
      f"{model_robust.translation[0]:.4f}), "
      f"Rotation: ({-model_robust.rotation:.4f}")


# visualise correspondences
fig, ax = plt.subplots(nrows=2, ncols=1)

plt.gray()

inliers_idxs = np.nonzero(inliers)[0]
plot_matches(ax[0], img_orig_gray, img_warped_gray, src, dst,
             np.column_stack((inliers_idxs, inliers_idxs)), matches_color='b')
ax[0].axis('off')
ax[0].set_title('Correct correspondences')

outlier_idxs = np.nonzero(outliers)[0]
plot_matches(ax[1], img_orig_gray, img_warped_gray, src, dst,
             np.column_stack((outlier_idxs, outlier_idxs)), matches_color='r')

ax[1].axis('off')
ax[1].set_title('Faulty correspondences')

plt.show()


#================================================================================
#                           Test on two Sentinel-2 scenes
#================================================================================

# read in Sentinel-2 GeoTiff
# reference image is affected by heavy cloud cover i.e. > 60%
# Can image co-registration still work under these conditions?
img_ref_fn = '/media/ubuntu/data_archive/F2020/Colombia/raster/cartagena/images/' \
        'S2A_MSIL2A_20180119T152631_N0206_R025_T18NWF_20180119T184927.tif'
# image to be co-registered almost cloud free
img_coreg_fn = '/media/ubuntu/data_archive/F2020/Colombia/raster/cartagena/images/' \
            'S2A_MSIL2A_20180208T152641_N0206_R025_T18NWF_20180208T185301.tif'
# read image in as array
img_ref = pyeo.raster_manipulation.raster_to_array(img_ref_fn)
img_coreg = pyeo.raster_manipulation.raster_to_array(img_coreg_fn)




def coregister_image():
    """
    This function co-registers a pair of images using scikit-image functions using RANSAC.

    Example taken from: https://scikit-image.org/docs/dev/auto_examples/transform/plot_matching.html

    :return: Path to co-registered image and image transform parameters.
    """




if __name__ == "__main__":




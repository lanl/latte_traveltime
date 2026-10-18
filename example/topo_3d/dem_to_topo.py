#
# Convert the DEM into the topography of the model
#
# The script reprojects the DEM to UTM zone 10N, smooths it, and samples it at the grid
# points of the model.
#
import numpy as np
import rasterio
from rasterio.transform import from_origin
from rasterio.warp import reproject, Resampling
from scipy.ndimage import gaussian_filter, map_coordinates
from pathlib import Path

# Model grid: 201 x 201 points with a spacing of 60 m along x (east) and y (north);
# the origin, x = y = 0, lies at 588500 m E and 5183500 m N in UTM zone 10N
# (EPSG:32610), so that the summit lies near the center of the model
nx = ny = 201
d = 60.0
x0 = 588500.0
y0 = 5183500.0

# Reproject the DEM to UTM on a 20 m grid that covers the model with a margin of 1960 m;
# the grid starts at its upper-left corner (west, north)
res = 20.0
margin = 1960.0
west, north = x0 - margin, y0 + (ny - 1)*d + margin
width = int(round(((nx - 1)*d + 2*margin)/res))
height = int(round(((ny - 1)*d + 2*margin)/res))
transform = from_origin(west, north, res, res)
dem = np.zeros((height, width), dtype=np.float32)
with rasterio.open('dem/rainier.tif') as src:
    reproject(source=rasterio.band(src, 1), destination=dem, src_transform=src.transform, src_crs=src.crs,
              dst_transform=transform, dst_crs='EPSG:32610', resampling=Resampling.cubic)
assert dem.min() > 0, 'The DEM does not cover the model'

# Smooth with a Gaussian of 120 m, twice the grid spacing of the model, to avoid
# aliasing when sampling at the grid points
dem = gaussian_filter(dem.astype(np.float64), 120.0/res)

# Sample at the grid points of the model; the center of pixel (row, column) lies at
# x = west + (column + 0.5)*res and y = north - (row + 0.5)*res
X, Y = np.meshgrid(np.arange(nx)*d, np.arange(ny)*d)
h = map_coordinates(dem, [(north - y0 - Y)/res - 0.5, (x0 + X - west)/res - 0.5], order=1)

s = np.degrees(np.arctan(np.hypot(*np.gradient(h, d))))
i, j = np.unravel_index(np.argmax(h), h.shape)
print(f'Elevation: {h.min():.0f} to {h.max():.0f} m, with the summit at x = {j*d:.0f} m, y = {i*d:.0f} m')
print(f'Slope: mean {s.mean():.1f} deg, 95th percentile {np.percentile(s, 95):.1f} deg, maximum {s.max():.1f} deg')

# Topography (elevation above sea level) as an (ny, nx) float32 array with y as the fastest axis
Path('model').mkdir(exist_ok=True)
h.T.astype(np.float32).tofile('model/topo.bin')

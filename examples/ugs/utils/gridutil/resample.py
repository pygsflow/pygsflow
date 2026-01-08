import rasterio
from rasterio import features
from rasterio.mask import mask
from shapely.geometry import Polygon
import numpy as np
from typing import List, Callable, Any
from rasterio.errors import RasterioIOError
import unittest
import tempfile
import os

def resample_raster_with_polygons(
    raster_path: str,
    polygons: List[Polygon],
    reduction_function: Callable[[np.ndarray], Any],
    nodata_value: float = np.nan,
) -> List[Any]:
    """Sample raster for each polygon using a reduction function."""
    sampled_values = []

    if reduction_function is None or reduction_function == "mean":
        reduction_function = np.mean   
    elif reduction_function == "median":
        reduction_function = np.median
    elif reduction_function == "max":
        reduction_function = np.max
    elif reduction_function == "min":
        reduction_function = np.min
    elif reduction_function == "sum":
        reduction_function = np.sum
    elif reduction_function == "std":
        reduction_function = np.std
    elif reduction_function == "var":
        reduction_function = np.var
    elif reduction_function == "count":
        reduction_function = np.size
    else:
       pass # use the reduction function as is

    try:
        with rasterio.open(raster_path) as src:
            for polygon in polygons:
                try:
                    out_image, out_transform = mask(src, [polygon], crop=True)
                except Exception as e:
                    sampled_values.append(np.nan)                        
                    continue
                
                nodata_value = src.nodata
                valid_pixels = out_image[out_image != nodata_value]

                if valid_pixels.size > 0:
                    # remove nodata values
                    valid_pixels = valid_pixels[valid_pixels != nodata_value]
                    result = reduction_function(valid_pixels)
                    sampled_values.append(result)
                else:
                    sampled_values.append(np.nan)
              
    except FileNotFoundError:
        print(f"Error: The raster file at '{raster_path}' was not found.")
        return []
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return []

    return sampled_values


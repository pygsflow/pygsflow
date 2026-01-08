import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from shapely.geometry import Polygon, Point
from shapely import vectorized
from tqdm import tqdm
import flopy

def compute_grid_intersections_rasterized(source_grid, target_grid, 
                                        resolution=100):
    """
    Vectorized version - much faster for large grids.
    Rasterize source and target grids into a common pixel grid, 
    then compute overlaps as a DataFrame (sparse format).
    
    Parameters
    ----------
    source_grid : list of shapely.Polygon
        Source grid cells.
    target_grid : list of shapely.Polygon
        Target grid cells.
    resolution : int
        Number of pixels along x/y per source cell.
    
    Returns
    -------
    df : pandas.DataFrame
        Columns: source, target, area
    """
    from shapely.vectorized import contains
    
    # Convert to numpy arrays for vectorization
    source_bounds = np.array([poly.bounds for poly in source_grid])
    target_bounds = np.array([poly.bounds for poly in target_grid])

    intersects = (
        (source_bounds[:, 0][:, None] <= target_bounds[:, 2]) &
        (source_bounds[:, 2][:, None] >= target_bounds[:, 0]) &
        (source_bounds[:, 1][:, None] <= target_bounds[:, 3]) &
        (source_bounds[:, 3][:, None] >= target_bounds[:, 1])
    )
    
    # Only process intersecting pairs
    source_indices, target_indices = np.where(intersects)
    
    total_pairs = len(source_grid) * len(target_grid)
    print(f"Found {len(source_indices)} potential intersections out of "
          f"{total_pairs} total pairs")
    
    records = []
    
    # Process only intersecting pairs with progress bar
    for idx, (s_idx, t_idx) in enumerate(
        tqdm(zip(source_indices, target_indices), 
              desc="Computing intersections", 
              total=len(source_indices))):
        source_poly = source_grid[s_idx]
        target_poly = target_grid[t_idx]
        
        # bounding box for source
        minx, miny, maxx, maxy = source_poly.bounds
        xs = np.linspace(minx, maxx, resolution)
        ys = np.linspace(miny, maxy, resolution)
        xx, yy = np.meshgrid(xs, ys)

        # rasterize source cell
        mask_source = vectorized.contains(source_poly, xx, yy)
        dx = (maxx - minx) / (resolution - 1)
        dy = (maxy - miny) / (resolution - 1)
        pixel_area = dx * dy

        # rasterize target cell on same pixel grid
        mask_target = vectorized.contains(target_poly, xx, yy)
        overlap = np.sum(mask_source & mask_target)

        if overlap > 0:            
            records.append({
                "source": s_idx,
                "target": t_idx,
                "area": overlap * pixel_area
            })

    return pd.DataFrame(records)


def compute_grid_exact_intersections(source_grid, target_grid):
    """
    Compute exact intersection areas between two grids using direct geometric 
    operations. Much faster and more accurate than rasterization approach.
    
    Parameters
    ----------
    source_grid : list of shapely.Polygon
        Source grid cells.
    target_grid : list of shapely.Polygon
        Target grid cells.
    
    Returns
    -------
    df : pandas.DataFrame
        Columns: source, target, area
    """
    # Convert to numpy arrays for vectorization
    source_bounds = np.array([poly.bounds for poly in source_grid])
    target_bounds = np.array([poly.bounds for poly in target_grid])
    
    # Vectorized intersection test using bounding boxes
    # This eliminates most non-overlapping pairs quickly
    intersects = (
        (source_bounds[:, 0][:, None] <= target_bounds[:, 2]) &
        (source_bounds[:, 2][:, None] >= target_bounds[:, 0]) &
        (source_bounds[:, 1][:, None] <= target_bounds[:, 3]) &
        (source_bounds[:, 3][:, None] >= target_bounds[:, 1])
    )
    
    # Only process intersecting pairs
    source_indices, target_indices = np.where(intersects)
    
    total_pairs = len(source_grid) * len(target_grid)
    print(f"Found {len(source_indices)} potential intersections out of "
          f"{total_pairs} total pairs")
    
    records = []
    
    # Process only intersecting pairs with progress bar
    for idx, (s_idx, t_idx) in enumerate(
        tqdm(zip(source_indices, target_indices), 
              desc="Computing exact intersections", 
              total=len(source_indices))):
        source_poly = source_grid[s_idx]
        target_poly = target_grid[t_idx]
        
        # Compute exact intersection area using Shapely
        intersection_poly = source_poly.intersection(target_poly)
        if not intersection_poly.is_empty:
            records.append({
                "source": s_idx,
                "target": t_idx,
                "area": intersection_poly.area
            })

    return pd.DataFrame(records)

def grid_to_shapely_polygons(grid):

    if isinstance(grid, flopy.discretization.StructuredGrid):
        polygons = []
        for i in range(grid.nrow):
            for j in range(grid.ncol):
                # Get cell vertices using FloPy's built-in method
                vertices = grid.get_cell_vertices(i, j)
                # Convert to 2D coordinates (ignore z if present)
                coords_2d = [(v[0], v[1]) for v in vertices]
                polygons.append(Polygon(coords_2d))
        return polygons
    elif isinstance(grid, flopy.discretization.UnstructuredGrid):
        polygons = []
        ncpl = int(grid.ncpl) if hasattr(grid, 'ncpl') else len(grid.iverts)
        for i in range(ncpl):
            # Get cell vertices for unstructured grid
            vertices = grid.get_cell_vertices(i)
            coords_2d = [(v[0], v[1]) for v in vertices]
            polygons.append(Polygon(coords_2d))
        return polygons
    elif isinstance(grid, flopy.discretization.VertexGrid):
        polygons = []        
        for i in range(len(grid.cell2d)):
            # Each row: [cell_num, x_center, y_center, n_vertices, v1, v2, v3, v4]
            cell_data = grid.cell2d[i]          
            n_vertices = int(cell_data[3])  # number of vertices
            
            # Extract vertex indices (columns 4 onwards)
            vertex_indices = [int(idx) for idx in cell_data[4:4+n_vertices] if int(idx) >= 0]
            
            # Get vertex coordinates from the vertices array            
            vertex_coords = []
            for idx in vertex_indices:
                vertex_coords.append((grid._vertices[idx][0], grid._vertices[idx][1]))
            
            polygons.append(Polygon(vertex_coords))
      
        return polygons
    else:
        # Fallback for unknown grid types
        raise ValueError(f"Unsupported grid type: {type(grid)}")


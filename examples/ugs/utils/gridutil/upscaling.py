import matplotlib.pyplot as plt
import numpy as np
from shapely.geometry import Polygon


def find_2x2_neighbors(index, rows, cols):
    """
    Finds the indices of a 2x2 block starting at the given index.
    Assumes a row-major grid.
    
    Args:
        index (int): The index of the top-left cell of the potential block.
        rows (int): The total number of rows in the grid.
        cols (int): The total number of columns in the grid.
        
    Returns:
        list: A list of the 3 other indices in the 2x2 block.
    """
    # Get the row and column from the index
    r = index // cols
    c = index % cols

    neighbors = []
    # Add the top-right neighbor
    if c + 1 < cols:
        neighbors.append(r * cols + (c + 1))
    # Add the bottom-left neighbor
    if r + 1 < rows:
        neighbors.append((r + 1) * cols + c)
    # Add the bottom-right neighbor
    if r + 1 < rows and c + 1 < cols:
        neighbors.append((r + 1) * cols + (c + 1))
        
    return neighbors

def coarsen_grid(fine_grid, p, rows, cols):
    """
    Coarsens a grid of polygons, keeping areas that intersect
    with polygon p at fine resolution.
    
    Args:
        fine_grid (list): The list of fine-resolution polygons.
        p (Polygon): The polygon defining the fine-resolution area.
        rows (int): The number of rows in the grid.
        cols (int): The number of columns in the grid.
        
    Returns:
        tuple: (coarsened_grid, mapping_table) where:
            - coarsened_grid: The new coarsened grid
            - mapping_table: List of tuples (fine_id, coarse_id)
    """
    coarsened_grid = []
    processed_indices = set()
    mapping_table = []
    
    # Process in 2x2 blocks to ensure proper coarsening
    for r in range(0, rows, 2):
        for c in range(0, cols, 2):
            # Get the 4 indices for the 2x2 block
            block_indices = []
            for dr in range(2):
                for dc in range(2):
                    if r + dr < rows and c + dc < cols:
                        block_indices.append((r + dr) * cols + (c + dc))
            
            # Skip if any cell in the block is already processed
            if any(idx in processed_indices for idx in block_indices):
                continue
            
            # Check if any cell in the block intersects with p
            intersects_p = any(fine_grid[idx].intersects(p) for idx in block_indices)
            
            if intersects_p:
                # Keep all cells in the block at fine resolution
                for idx in block_indices:
                    coarse_id = len(coarsened_grid)
                    coarsened_grid.append(fine_grid[idx])
                    mapping_table.append((idx, coarse_id))
                    processed_indices.add(idx)
            else:
                # All cells are outside p, can be coarsened
                if len(block_indices) == 4:
                    # Full 2x2 block - merge into one polygon
                    polygons_to_merge = [fine_grid[idx] for idx in block_indices]
                    coarsened_polygon = polygons_to_merge[0]
                    for poly in polygons_to_merge[1:]:
                        coarsened_polygon = coarsened_polygon.union(poly)
                    
                    coarse_id = len(coarsened_grid)
                    coarsened_grid.append(coarsened_polygon)
                    
                    # Map all 4 fine cells to the same coarse cell
                    for idx in block_indices:
                        mapping_table.append((idx, coarse_id))
                else:
                    # Partial block (edge case) - keep individual cells
                    for idx in block_indices:
                        coarse_id = len(coarsened_grid)
                        coarsened_grid.append(fine_grid[idx])
                        mapping_table.append((idx, coarse_id))
                
                # Mark all cells in the block as processed
                for idx in block_indices:
                    processed_indices.add(idx)
    
    # Handle any remaining unprocessed cells (shouldn't happen with proper 2x2 processing)
    for i in range(len(fine_grid)):
        if i not in processed_indices:
            coarse_id = len(coarsened_grid)
            coarsened_grid.append(fine_grid[i])
            mapping_table.append((i, coarse_id))
            processed_indices.add(i)

    return coarsened_grid, mapping_table

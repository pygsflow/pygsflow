import matplotlib.pyplot as plt
import numpy as np
from shapely.geometry import Polygon
from typing import List, Optional, Union
from .grid import Grid

def plot_grid(polygons: List[Polygon], 
              values: Optional[List[float]] = None,
              title: str = "Grid Plot",
              figsize: tuple = (10, 8),
              ax: Optional[plt.Axes] = None,
              show_ids: bool = False,
              colormap: str = 'viridis',
              alpha: float = 0.7,
              edge_color: str = 'black',
              edge_width: float = 0.5,
              show_colorbar: bool = True) -> None:
    """
    Plot a grid represented as a list of Shapely polygons.
    
    Args:
        polygons (List[Polygon]): List of Shapely polygon objects 
                                 representing grid cells
        values (Optional[List[float]]): Optional values to color the 
                                       polygons by. If None, polygons 
                                       are colored by index.
        title (str): Title for the plot
        figsize (tuple): Figure size (width, height)
        show_ids (bool): Whether to show cell IDs as text in the center 
                        of each polygon
        colormap (str): Matplotlib colormap name
        alpha (float): Transparency of polygon fills (0-1)
        edge_color (str): Color of polygon edges
        edge_width (float): Width of polygon edges
        show_colorbar (bool): Whether to show colorbar when values are 
                             provided
    """
    # if ax is passed, use it, otherwise create a new figure
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()
    
    # Determine colors
    if values is not None:
        if len(values) != len(polygons):
            raise ValueError(f"Number of values ({len(values)}) must "
                           f"match number of polygons ({len(polygons)})")
        
        # Normalize values for colormap
        values = np.array(values)
        if np.any(np.isnan(values)):
            # Handle NaN values
            valid_mask = ~np.isnan(values)
            if not np.any(valid_mask):
                # All values are NaN, use uniform color
                colors = ['lightgray'] * len(polygons)
            else:
                # Create colors with NaN handling
                valid_values = values[valid_mask]
                cmap = plt.cm.get_cmap(colormap)
                # Normalize valid values
                normalized_values = ((valid_values - valid_values.min()) / 
                                   (valid_values.max() - valid_values.min()))
                # Create color array
                colors = ['lightgray'] * len(polygons)
                valid_idx = 0
                for i, is_valid in enumerate(valid_mask):
                    if is_valid:
                        colors[i] = cmap(normalized_values[valid_idx])
                        valid_idx += 1
        else:
            # No NaN values, normal coloring
            cmap = plt.cm.get_cmap(colormap)
            normalized_values = ((values - values.min()) / 
                               (values.max() - values.min()))
            colors = [cmap(v) for v in normalized_values]
    else:
        # # Color by index
        # cmap = plt.cm.get_cmap(colormap)
        # colors = [cmap(i / len(polygons)) for i in range(len(polygons))]
        colors = None
    
    # Plot polygons
    for i, polygon in enumerate(polygons):
        if colors is not None:            
            x, y = polygon.exterior.xy
            ax.fill(x, y, alpha=alpha, edgecolor=edge_color, 
                    linewidth=edge_width, facecolor=colors[i])
        else:
            x, y = polygon.exterior.xy
            ax.fill(x, y, alpha=alpha, edgecolor=edge_color, 
                    linewidth=edge_width, facecolor=None)

        
        # Add cell ID text if requested
        if show_ids:
            centroid = polygon.centroid
            ax.text(centroid.x, centroid.y, str(i), ha='center', 
                   va='center', fontsize=8, fontweight='bold', 
                   color='white' if values is None else 'black')
    
    # Add colorbar if values are provided and requested
    if (values is not None and show_colorbar and 
        not np.all(np.isnan(values))):
        valid_values = values[~np.isnan(values)]
        if len(valid_values) > 0:
            sm = plt.cm.ScalarMappable(cmap=colormap, 
                                     norm=plt.Normalize(
                                         vmin=valid_values.min(), 
                                         vmax=valid_values.max()))
            sm.set_array([])
            cbar = plt.colorbar(sm, ax=ax)
            cbar.set_label('Values')
    
    ax.set_title(title)
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_aspect('equal')    
    plt.tight_layout()
    

def _plot_grid_on_axis(ax, polygons: List[Polygon], 
                      values: Optional[List[float]] = None,
                      title: str = "Grid",
                      show_ids: bool = False,
                      colormap: str = 'viridis') -> None:
    """Helper function to plot grid on a specific axis."""
    # Determine colors
    if values is not None:
        if len(values) != len(polygons):
            raise ValueError(f"Number of values ({len(values)}) must "
                           f"match number of polygons ({len(polygons)})")
        
        values = np.array(values)
        if np.any(np.isnan(values)):
            valid_mask = ~np.isnan(values)
            if not np.any(valid_mask):
                colors = ['lightgray'] * len(polygons)
            else:
                valid_values = values[valid_mask]
                cmap = plt.cm.get_cmap(colormap)
                normalized_values = ((valid_values - valid_values.min()) / 
                                   (valid_values.max() - valid_values.min()))
                colors = ['lightgray'] * len(polygons)
                valid_idx = 0
                for i, is_valid in enumerate(valid_mask):
                    if is_valid:
                        colors[i] = cmap(normalized_values[valid_idx])
                        valid_idx += 1
        else:
            cmap = plt.cm.get_cmap(colormap)
            colors = (values - values.min()) / (values.max() - values.min())
            colors = [cmap(c) for c in colors]
    else:
        cmap = plt.cm.get_cmap(colormap)
        colors = [cmap(i / len(polygons)) for i in range(len(polygons))]
    
    # Plot polygons
    for i, polygon in enumerate(polygons):
        x, y = polygon.exterior.xy
        ax.fill(x, y, alpha=0.7, edgecolor='black', linewidth=0.5, 
                facecolor=colors[i])
        
        if show_ids:
            centroid = polygon.centroid
            ax.text(centroid.x, centroid.y, str(i), ha='center', 
                   va='center', fontsize=8, fontweight='bold', 
                   color='white')
    
    ax.set_title(title)
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_aspect('equal')

def plot_layer(grid: Grid, layer: int, title: str = "Layer Plot", 
               figsize: tuple = (10, 8), ax: Optional[plt.Axes] = None, 
               show_ids: bool = False, colormap: str = 'viridis', 
               alpha: float = 0.7, edge_color: str = 'black', 
               edge_width: float = 0.5, show_colorbar: bool = True) -> None:
    """Plot a layer of a grid."""
    # from grid get geometry, which is a list of polygons.
    polygons = grid.geometry
    polys = polygons['geometry'].tolist()
    plot_grid(polys)
    cc = 1
    pass

def fast_plotter_3d(grid: Grid, layer: int, title: str = "Layer Plot") -> None:
    import pyvista as pv
    import numpy as np
    
    # Get points as numpy array
    points = np.array(grid.vertices)
    
    # Get vertex indices for the specified layer
    layer_mask = grid.node_data['lay'] == layer
    layer_vertex_indices = grid.node_data[layer_mask]['vertex_indices']
    
    # Convert to proper format for PyVista
    # PyVista expects separate arrays for cells and cell types
    cells = []
    cell_types = []
    
    for vertex_list in layer_vertex_indices:
        # Each cell: [n_vertices, v1, v2, ..., v8]
        cells.extend([8] + [int(v-1) for v in vertex_list])
        # Cell type 12 is VTK_HEXAHEDRON
        cell_types.append(12)
    
    cells = np.array(cells, dtype=np.int32)
    cell_types = np.array(cell_types, dtype=np.uint8)
    
    # Create unstructured grid
    ugrid = pv.UnstructuredGrid(cells, cell_types, points)
    
    # Create plotter
    plotter = pv.Plotter()
    plotter.add_mesh(ugrid, show_edges=True, cmap="viridis")
    plotter.view_xy()  # top-down 2D view
    plotter.show(title=f"Layer {layer} - {title}")
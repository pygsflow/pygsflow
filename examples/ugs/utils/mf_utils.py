import numpy as np
import pandas as pd
from shapely.geometry import LineString, Point
from shapely.geometry import Polygon, LineString
from shapely.strtree import STRtree
from gridutil import plot_grid
import matplotlib.pyplot as plt


def get_sfr_lines(mf):
    reach_data = mf.sfr.reach_data
    df_reach_data = pd.DataFrame(reach_data)
    
    cell_size = mf.modelgrid.delr[0]
    lines = []
    line_info = []
    for segment in df_reach_data['iseg'].unique():
        df_segment = df_reach_data[
            df_reach_data['iseg'] == segment
        ].sort_values('ireach')
        
        # Create lines connecting consecutive reaches
        #seg_line =[]
        for i in range(len(df_segment) - 1):
            current_reach = df_segment.iloc[i]
            next_reach = df_segment.iloc[i + 1]
            
            # Calculate center coordinates
            current_center = (
                (current_reach['j']-1) * cell_size + cell_size / 2,
                (current_reach['i']-1) * cell_size + cell_size / 2
            )
            next_center = (
                (next_reach['j']-1) * cell_size + cell_size / 2,
                (next_reach['i']-1) * cell_size + cell_size / 2
            )
            
            line = LineString([current_center, next_center])
            #lines.append(line)
            line_info.append([segment, current_reach['ireach'], next_reach['ireach'], line])
    lines_df = pd.DataFrame(line_info, columns=['segment', 'ireach', 'next_ireach', 'lineGeometry'])
    return lines_df


def get_unique_vertices(dis):
    """
    Find unique vertices in a structured grid.
    
    Parameters:
    -----------
    dis : flopy discretization object
        Discretization package containing grid information
    
    Returns:
    --------
    vertices : numpy.ndarray
        Array of unique vertex coordinates (x, y, z)
    """
    delr = dis.delr.array
    delc = dis.delc.array
    nrow = dis.nrow
    ncol = dis.ncol
    nlay = dis.nlay
    top = dis.top.array
    bottom = dis.botm.array
    
    # Pre-calculate coordinate arrays
    x_coords = np.concatenate([[0], np.cumsum(delr)])
    y_coords = np.concatenate([[0], np.cumsum(delc)])
    
    # Pre-allocate arrays for better performance
    max_vertices = (nrow + 1) * (ncol + 1) * (nlay + 2)
    vertices_x = np.empty(max_vertices, dtype=np.float64)
    vertices_y = np.empty(max_vertices, dtype=np.float64)
    vertices_z = np.empty(max_vertices, dtype=np.float64)
    
    vertex_count = 0
    
    # Vectorized approach for better performance
    for i in range(nrow + 1):
        for j in range(ncol + 1):
            x = x_coords[j]
            y = y_coords[i]
            
            # Collect z values for this location
            z_values = []
            
            # Add top elevation if within grid bounds
            if i < nrow and j < ncol:
                z_values.append(top[i, j])
            
            # Add bottom elevations if within grid bounds
            if i < nrow and j < ncol:
                z_values.extend(bottom[:, i, j])
            
            # Remove duplicates and sort
            z_values = np.unique(z_values)
            
            # Add vertices for this (x,y) location
            n_z = len(z_values)
            vertices_x[vertex_count:vertex_count + n_z] = x
            vertices_y[vertex_count:vertex_count + n_z] = y
            vertices_z[vertex_count:vertex_count + n_z] = z_values
            vertex_count += n_z
    
    # Trim arrays to actual size and combine
    vertices = np.column_stack([
        vertices_x[:vertex_count],
        vertices_y[:vertex_count],
        vertices_z[:vertex_count]
    ])
    
    return vertices

def intersect_line_with_grid(mi):    
    lines = mi.sfr_df["seg_line"].values
    polygons = mi.oct_grid2d

    tree = STRtree(polygons)

    usg_sfr = []
    sfr_polygons = []
   
    # 2. Query the tree for each line
    for i, line in enumerate(lines):
        # tree.query(line) efficiently returns *only* the geometry objects 
        # from 'polygons' that intersect 'line'
        intersecting_geometries = tree.query(line, predicate='intersects')

        # Create a list of the *indices* of the intersecting polygons
        intersecting_indices = []

        # Map the resulting geometries back to their original index
        # Note: This is where floating point precision could still be an issue 
        # if the query results are not exactly equal to the polygons list members.
        min_length_threshold = 5
        x0 = line.xy[0][0]
        y0 = line.xy[1][0]
        segment_start_point = Point(x0, y0)
        debug = False
        if debug:
            x, y = line.xy
            plt.plot(x, y, 'k-')
            plt.show()
            pass

        for poly in intersecting_geometries:          
                     
            intersection_geom = line.intersection(polygons[poly])
            if intersection_geom.length >= min_length_threshold:
                sfr_polygons.append(polygons[poly])
                intersecting_indices.append(f"Polygon {poly}")
                p = polygons[poly]
                p_centroid = p.centroid
                distance = segment_start_point.distance(p_centroid)
                # segment_id, polygon_id, distance form segment_upstream, intersection_length                   
                seg_inter_ = [i+1, poly, distance, intersection_geom.length ]
                usg_sfr.append(seg_inter_)
                
            

    mi.usg_sfr_df = pd.DataFrame(usg_sfr, columns=['segment_id', 'node_id', 'distance_from_segment_upstream', 'intersection_length'])    
    Debug = False
    if Debug:
        # plots the lines and the polygons
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()       
        
        plot_grid(polygons, ax=ax, edge_color='green')
        plot_grid(sfr_polygons, ax=ax, edge_color='red')
        for line in lines:
            ax.plot(line.xy[0], line.xy[1], 'k-')
       
        plt.show()
   
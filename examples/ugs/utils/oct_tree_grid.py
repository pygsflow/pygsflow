import os
import sys
from pprint import pformat
from tempfile import TemporaryDirectory
from itertools import product

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from shapely.geometry import Polygon, LineString

# run installed version of flopy or add local path
import flopy
from flopy.utils import flopy_io
from flopy.utils.gridgen import Gridgen
from utils import prms_utils

from utils import gridutil


# ================================
# (1) Assert Gridgen Executable
# ================================
def assert_gridgen_executable():
    gridgen_exe = flopy.which("gridgen")
    if gridgen_exe is None:
        msg = (
            "Warning, gridgen is not in your path. "
            "When you create the griden object you will need to "
            "provide a full path to the gridgen binary executable."
        )
        print(msg)
    else:
        print(f"gridgen executable was found at: "
              f"{flopy_io.relpath_safe(gridgen_exe)}")

def get_line_coords(sfr):
    """
    Extract line coordinates from a stream object for Gridgen refinement.
    Uses stream topology (outseg) to properly connect segments.
    
    Parameters
    ----------
    stream_obj : _StreamsObj
        Stream object from FlowAccumulation.make_streams()
    modelgrid : StructuredGrid
        Model grid object with coordinate information
    
    Returns
    -------
    tuple
        (stream_lines, ij_coords) where stream_lines is list of line 
        coordinates and ij_coords is list of ij coordinates for Gridgen
    """
    # Get stream cells (where iseg > 0 indicates a stream cell)
    nrow = sfr.parent.nrow
    ncol = sfr.parent.ncol
    cell_size = sfr.parent.modelgrid.delr[0]
    delr = sfr.parent.modelgrid.delr
    delc = sfr.parent.modelgrid.delc
    segment_data = pd.DataFrame(sfr.segment_data[0])
    reach_data = pd.DataFrame(sfr.reach_data)
    modelgrid = sfr.parent.modelgrid
    unique_iseg = segment_data["nseg"].unique()

    #x_coords, y_coords = modelgrid.xycenters
  
    reach_info = []
    for iseg in unique_iseg:               
        if iseg > 0:
            curr_seg_info =  reach_data[reach_data["iseg"] == iseg]
            reach_ids = curr_seg_info["ireach"].values
            seg_line = []
            ij_coords = []
            reach_ids = np.sort(reach_ids)
            for ireach in reach_ids:               
                i = curr_seg_info[curr_seg_info["ireach"] == ireach]["i"].values
                j = curr_seg_info[curr_seg_info["ireach"] == ireach]["j"].values              
                x1 = delr[j]*j - delr[j]/2
                y1 = delc.sum() - (delc[i]*i - delc[i]/2)
                eps = cell_size/10000          
                seg_line.append((x1+eps, y1+eps))
                ij_coords.append((i[0], j[0]))
            outseg = segment_data[segment_data["nseg"] == iseg]["outseg"].values[0]            
            reach_info.append([iseg, seg_line, outseg, ij_coords])
        
    df = pd.DataFrame(reach_info, columns=["iseg", "seg_line", "outseg", "ij_coords"])
    # connect segments
    for iseg in unique_iseg:
        df_seg = df[df["iseg"] == iseg]
        df_seg_line = df_seg["seg_line"].values[0]
        outseg = df_seg["outseg"].values[0]

        if outseg > 0:
            df_outseg = df[df["iseg"] == outseg]
            df_outseg_line = df_outseg["seg_line"].values[0]
            first_point = df_outseg_line[0]
            df_seg_line.append(first_point)

        df.loc[df["iseg"] == iseg, "seg_line"] = LineString(df_seg_line)
        
    
    debug = False
    if debug:
        fig, ax = plt.subplots()
        for i, row in df.iterrows():
            ax.plot(row["seg_line"].xy[0], row["seg_line"].xy[1], 'k-')
        plt.show()
    return df



   
      

def mask_to_polygon(mask, modelgrid):
    """
    Convert a boolean mask to a polygon boundary using matplotlib contours.
    Returns the polygon in Gridgen format: [[[vertices]]]
    
    Parameters
    ----------
    mask : numpy.ndarray
        Boolean mask array (True for watershed, False for outside)
    modelgrid : StructuredGrid
        Model grid object with coordinate information
    
    Returns
    -------
    list
        Polygon in Gridgen format: [[[vertices]]] - feature collection 
        containing polygon list containing vertex list
    """
    # Create a figure and axis for contour plotting
    fig, ax = plt.subplots()
    
    # Use row and column numbers instead of coordinates
    nrow, ncol = mask.shape
    row_numbers = np.arange(nrow)
    col_numbers = np.arange(ncol)
    
    # Create meshgrid for row and column indices
    row_grid, col_grid = np.meshgrid(row_numbers, col_numbers, indexing='ij')
    row_grid = row_grid * modelgrid.delr[0]
    col_grid = col_grid * modelgrid.delc[0]
    
    # Find contours of the mask (boundary between True and False)
    contours = ax.contour(col_grid, np.flipud(row_grid), mask.astype(float), levels=[0.5])
    
    # Get the contour paths - QuadContourSet has allsegs attribute
    if hasattr(contours, 'allsegs') and contours.allsegs:
        # Get the largest contour (main watershed boundary)
        all_paths = contours.allsegs[0]  # First level (0.5)
        if all_paths:
            largest_path = max(all_paths, key=lambda p: len(p))
            # Convert each vertex from [x, y] to (x, y) tuple format
            polygon_coords = [tuple(vertex) for vertex in largest_path.tolist()]
        else:
            polygon_coords = []
    else:
        polygon_coords = []
    
    # Close the figure to free memory
    plt.close(fig)
    
    if polygon_coords:
        return [[polygon_coords]]
    else:
        print("Warning: No contour found for watershed mask")
        return []

def create_buffer_polygons_around_streams(stream_lines, buffer_distance, 
                                          cell_size, 
                                          active_domain_polygon=None, 
                                          inset_distance=10.0):
    """
    Create a single dissolved buffer polygon around all stream lines for Gridgen refinement.
    The buffer is clipped to stay within the active domain with an inset from the boundary.
    
    Parameters
    ----------
    stream_lines : list
        List of stream line coordinates
    buffer_distance : float
        Buffer distance in model units
    cell_size : float
        Size of grid cells for scaling
    active_domain_polygon : list, optional
        Active domain polygon in Gridgen format to clip the buffer against
    inset_distance : float, optional
        Distance to inset from the active domain boundary (default: 10.0)
    
    Returns
    -------
    list
        List containing a single dissolved buffer polygon in Gridgen format
    """
    try:
        from shapely.geometry import LineString, Polygon
        from shapely.ops import unary_union
        
        individual_polygons = []
        
        for line in stream_lines:         
            buffered = line.buffer(buffer_distance)
            individual_polygons.append(buffered)
        
        if not individual_polygons:
            return []
        
        # Dissolve all polygons into one using unary_union
        dissolved_polygon = unary_union(individual_polygons)
        
        # Clip buffer polygon against active domain if provided
        if active_domain_polygon and len(active_domain_polygon) > 0:
            # Convert active domain polygon to Shapely polygon
            active_domain_coords = active_domain_polygon[0][0]  # Get first polygon coordinates
            active_domain_shapely = Polygon(active_domain_coords)
            
            # Create inset polygon by shrinking the active domain
            inset_polygon = active_domain_shapely.buffer(-inset_distance)
            
            # Intersect buffer with inset active domain
            if not inset_polygon.is_empty:
                clipped_polygon = dissolved_polygon.intersection(inset_polygon)
                
                # Use clipped polygon if it's not empty
                if not clipped_polygon.is_empty:
                    dissolved_polygon = clipped_polygon
        
        # Convert dissolved polygon to Gridgen format
        if hasattr(dissolved_polygon, 'exterior'):
            # Single polygon
            coords = list(dissolved_polygon.exterior.coords)
            return [[coords]]
        else:
            # MultiPolygon - extract each polygon
            dissolved_polygons = []
            for geom in dissolved_polygon.geoms:
                coords = list(geom.exterior.coords)
                dissolved_polygons.append([coords])
            return dissolved_polygons
        
    except ImportError:
        print("Warning: Shapely not available, using simple rectangular "
              "buffers")
        return create_simple_buffer_polygons(stream_lines, buffer_distance, 
                                            cell_size, 
                                            active_domain_polygon, 
                                            inset_distance)

def create_simple_buffer_polygons(stream_lines, buffer_distance, cell_size, 
                                  active_domain_polygon=None, 
                                  inset_distance=10.0):
    """
    Create a single dissolved rectangular buffer polygon around all stream 
    lines without Shapely.
    The buffer is clipped to stay within the active domain with an inset from the boundary.
    
    Parameters
    ----------
    stream_lines : list
        List of stream line coordinates
    buffer_distance : float
        Buffer distance in model units
    cell_size : float
        Size of grid cells for scaling
    active_domain_polygon : list, optional
        Active domain polygon in Gridgen format to clip the buffer against
    inset_distance : float, optional
        Distance to inset from the active domain boundary (default: 10.0)
    
    Returns
    -------
    list
        List containing a single dissolved buffer polygon in Gridgen format
    """
    if not stream_lines:
        return []
    
    buffer_cells = buffer_distance / cell_size
    
    # Find the overall bounding box of all streams
    all_x_coords = []
    all_y_coords = []
    
    for line_coords in stream_lines:
        if len(line_coords) < 2:
            continue
        all_x_coords.extend([coord[0] for coord in line_coords])
        all_y_coords.extend([coord[1] for coord in line_coords])
    
    if not all_x_coords or not all_y_coords:
        return []
    
    # Create a single rectangular buffer around all streams
    min_x = min(all_x_coords) - buffer_cells
    max_x = max(all_x_coords) + buffer_cells
    min_y = min(all_y_coords) - buffer_cells
    max_y = max(all_y_coords) + buffer_cells
    
    # Create single rectangular polygon
    rect_coords = [
        (min_x, min_y),
        (max_x, min_y),
        (max_x, max_y),
        (min_x, max_y),
        (min_x, min_y)  # Close the polygon
    ]
    
    # If active domain is provided, clip the rectangular buffer
    if active_domain_polygon and len(active_domain_polygon) > 0:
        try:
            from shapely.geometry import Polygon
            
            # Convert active domain polygon to Shapely polygon
            active_domain_coords = active_domain_polygon[0][0]  # Get first polygon coordinates
            active_domain_shapely = Polygon(active_domain_coords)
            
            # Create inset polygon by shrinking the active domain
            inset_polygon = active_domain_shapely.buffer(-inset_distance)
            
            # Convert buffer rectangle to Shapely polygon
            buffer_shapely = Polygon(rect_coords)
            
            # Intersect buffer with inset active domain
            if not inset_polygon.is_empty:
                clipped_polygon = buffer_shapely.intersection(inset_polygon)
                
                # Use clipped polygon if it's not empty
                if not clipped_polygon.is_empty:
                    if hasattr(clipped_polygon, 'exterior'):
                        # Single polygon
                        coords = list(clipped_polygon.exterior.coords)
                        return [[coords]]
                    else:
                        # MultiPolygon - extract each polygon
                        clipped_polygons = []
                        for geom in clipped_polygon.geoms:
                            coords = list(geom.exterior.coords)
                            clipped_polygons.append([coords])
                        return clipped_polygons
        except ImportError:
            # If Shapely is not available, return the unclipped rectangle
            print("Warning: Cannot clip buffer without Shapely, returning "
                  "unclipped rectangle")
    
    return [[rect_coords]]

def grid_2d_from_gridgen(mi):
    g = mi.gridgen_obj
    mi.gridprops = g.get_gridprops_disu5() 
    mi.unstructured_grid = g.get_gridprops_unstructuredgrid()

    vertices = np.array(mi.unstructured_grid['vertices'])
    grid2d = []
    for icell, cell in enumerate(mi.unstructured_grid['iverts']):
        ver_idx = cell[:-1]
        poly = [(x,y) for x,y in vertices[ver_idx, 1:]]     
        poly = Polygon(poly)
        grid2d.append(poly)
    return grid2d
def grid_from_structured_dis(dis) :
    cols = ["node_id", "geometry", "top", "botm", "layer"]
    rows = []
    modelgrid = dis.parent.modelgrid
    for k in range(dis.nlay):
        for i in range(dis.nrow):
            for j in range(dis.ncol):
                node_id = modelgrid.get_node((k, i, j))
                vertices = modelgrid.get_cell_vertices(node_id)
                poly = Polygon(vertices)
                top = modelgrid.top[i, j]
                botm = modelgrid.botm[k, i, j]
                rows.append({
                    "node_id": node_id[0],
                    "geometry": poly,
                    "top": top,
                    "botm": botm,
                    "layer": k,
                })
    return pd.DataFrame(rows, columns=cols)
    #   grid = cls(
    #         delr=delr_array,
    #         delc=delc_array,
    #         nrows=dis.nrow,
    #         ncols=dis.ncol,
    #         nlays=dis.nlay,
    #         top=getattr(dis, 'top', 0.0),
    #         botm=getattr(dis, 'botm', 0.0),
    #         origin_x=getattr(dis, 'origin_x', 0.0),
    #         origin_y=getattr(dis, 'origin_y', 0.0),
    #         rotation=getattr(dis, 'rotation', 0.0)
    #     )
    #     grid.nlays = dis.nlay

def create_oct_tree_grid(mi, stream_buffer_distance=100.0):
    assert_gridgen_executable()   
    gridgen_ws = os.path.join(mi.ws, "gridgen")
    if not os.path.exists(gridgen_ws):
        os.makedirs(gridgen_ws, exist_ok=True)
    print(f"Model workspace is : {flopy_io.scrub_login(mi.ws)}")
    print(f"Gridgen workspace is : {flopy_io.scrub_login(gridgen_ws)}")

   
    ms = mi.fine_gsf.mf
    dis = mi.fine_gsf.mf.dis
    source_grid = grid_from_structured_dis(mi.fine_gsf.mf.dis)
    gridgen_exe = mi.gridgen_exe
    #ms.modelgrid.set_coord_info(crs=mi.fine_gsf.mf.crs)
    g = Gridgen(ms.modelgrid, model_ws=gridgen_ws, exe_name= gridgen_exe)

    # Convert watershed mask to polygon boundary
    mi.watershed = prms_utils.get_prms_watershed(mi.fine_gsf)
    polygon = mask_to_polygon(mi.watershed, mi.fine_gsf.mf.modelgrid)  

    g.add_active_domain(polygon, range(mi.nlay))

    # from the stream object, get line coordinates.
    sfr = mi.fine_gsf.mf.sfr
    sfr_df = get_line_coords(sfr)
    mi.sfr_df = sfr_df
    ij_coords = sfr_df["seg_line"].values.tolist()  
       
    # Create buffer polygons around streams for additional refinement
    cell_size = mi.fine_gsf.mf.modelgrid.delr[0]  # Assuming uniform cell size
    stream_buffer_distance= 60
    inset_distance = 10.0  # Distance to inset from watershed boundary
    # todo: when the polygon has too small external areas, it fails
    buffer_polygons = create_buffer_polygons_around_streams(ij_coords, 
                                                           stream_buffer_distance, 
                                                           cell_size, 
                                                           polygon, 
                                                           inset_distance)
    
    if buffer_polygons:
        print(f"Created {len(buffer_polygons)} buffer polygons around "
              f"streams")        
        g.add_refinement_features(buffer_polygons, "polygon", 2, 
                                  range(mi.nlay))   

    # Add stream lines with fine refinement, todo: issues in connect stream lines
    #g.add_refinement_features(ij_coords, "line", 3, range(mi.nlay))

    # build the grid
    g.build(verbose=True)  
    if 0:
        fig = plt.figure(figsize=(15, 15))
        ax = fig.add_subplot(1, 1, 1, aspect="equal")
        g.plot(ax, linewidth=0.5)


    mi.gridgen_obj = g
    mi.oct_grid2d = grid_2d_from_gridgen(mi)


    # create mapping array
    struct_grid = source_grid[source_grid['layer'] == 0]['geometry'].values.tolist()  

    mi.mapping_df = gridutil.compute_grid_exact_intersections(source_grid = struct_grid, target_grid = mi.oct_grid2d)
    
    # creat gsf file
    vertices = pd.DataFrame(np.array(mi.unstructured_grid['vertices']), columns=['id', 'x', 'y'])
    n_vertices_top = vertices.shape[0]
    node_data = pd.DataFrame(columns=['inode', 'x', 'y', 'z', 'lay', 'vertex_indices'])
    node_data['inode'] = 1+np.arange(mi.gridprops['nodes'])
    node_data['x'] = mi.unstructured_grid['xcenters']
    node_data['y'] = mi.unstructured_grid['ycenters']
    node_data['z'] = (mi.unstructured_grid['top'] + mi.unstructured_grid['botm']) / 2.0
    node_data['lay'] = 1
    vertices2 = []
    iverts2 = []
    vertex_counter = 1
    vertex_dict = {}
    if 0:
        for jj, iv in enumerate(mi.unstructured_grid['iverts']):
            df_ = vertices.iloc[iv[:-1], :]
            xy = df_[['x', 'y']].values.tolist()
            z1 = mi.unstructured_grid['top'][jj]
            z2 = mi.unstructured_grid['botm'][jj]
            part1 = []
            part2 = []
            for xy_ in xy:
                part1.append([xy_[0], xy_[1], z1])
                part2.append([xy_[0], xy_[1], z2])
            vertices2.extend(part1)
            vertices2.extend(part2)
            v_id = np.arange(vertex_counter, vertex_counter + 2*len(xy))
            v_id = v_id.tolist()
            iverts2.append(v_id)
            vertex_counter = vertex_counter + 2*len(xy)
    else:
       
        veritices_top = pd.DataFrame(np.array(mi.unstructured_grid['vertices']), columns=['id', 'x', 'y'])
        veritices_top['z'] = None
        veritices_top['id']=veritices_top['id'].astype(int)
        veritices_top.sort_values(by = 'id', inplace = True)

        veritices_bot = pd.DataFrame(np.array(mi.unstructured_grid['vertices']), columns=['id', 'x', 'y'])
        veritices_bot['z'] = None
        veritices_bot['id']=veritices_bot['id'].astype(int)
        veritices_bot.sort_values(by = 'id', inplace = True)

        for jj, iv in enumerate(mi.unstructured_grid['iverts']):
            df_ = vertices.iloc[iv[:-1], :]
            xy = df_[['x', 'y']].values.tolist()
            z1 = mi.unstructured_grid['top'][jj]
            z2 = mi.unstructured_grid['botm'][jj]

            veritices_top.loc[iv[:-1], 'z'] = z1
            veritices_bot.loc[iv[:-1], 'z'] = z2
            # ids_top = []
            # for xy_ in xy:
            #     mask = (veritices_top['x'] == xy_[0]) &  (veritices_top['y'] == xy_[1])
            #     veritices_top.loc[mask, 'z'] = z1
            #     veritices_bot.loc[mask, 'z'] = z2
            #     id_ = veritices_top.loc[mask, 'id'].values[0]
            #     ids_top.append(id_+1)

            ids_bot = [v + n_vertices_top for v in iv[:-1]]
            new_iverts = iv[:-1]+ids_bot
            new_iverts = [v+1 for v in new_iverts]
            # v_id_top = np.arange(vertex_counter, vertex_counter +len(xy))
            # v_id_bot = n_vertices_top + np.arange(vertex_counter, vertex_counter +len(xy))            
            iverts2.append(new_iverts)
            

    vertices2 = pd.concat([veritices_top, veritices_bot])
    node_data['vertex_indices'] = iverts2
    vertices2 = pd.DataFrame(vertices2, columns=['x', 'y', 'z'])   

    #node_data.loc[i, 'vertex_indices'] = i

   
    gridutil.write_gsf_file(vertices2, node_data, filename=r"C:\workspace\projects\gsflow6\usg\grid_spec.gsf")



    
    

    end = 1

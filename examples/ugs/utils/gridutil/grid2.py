from shapely.geometry import Polygon
import tempfile
import os
import numpy as np
import pandas as pd
import geopandas as gpd
from typing import List


class Grid():
    """
    Base class for grid objects.
    Assumptions:
    - vertices are in the order of the node_data dataframe    
    - The grid prism in 3D, meaning that for any cell top and bottom 
      faces are parallel and have the same 2D shape
    - the attribute geometry is a dataframe with columns inode, 2d 
      shapely polygons, z_top and z_bottom)

    """
    def __init__(self):  
       
        self._nnode = None
       
        self.vertices = []  # dataframe with columns x, y
        self.node_data = []  # dataframe with columns inode, x, y, z, lay, 
                             # vertex_indices
        self._geometry = None  # dataframe with columns node_id, geometry, 
                               # top, botom
        
        self._nlays = None
        
    @property
    def nlays(self):
        if not (self._nlays is None):
            return self._nlays

        cell_1 = self.node_data[self.node_data['inode'] == 0]
        x = cell_1['x'].values[0]
        y = cell_1['y'].values[0]
        cells = self.node_data[(self.node_data['x'] == x) & 
                               (self.node_data['y'] == y)]
        z = cells['z'].unique()
        self._nlays = z.shape[0]        
        return self._nlays
    
    @property
    def nnode(self):
        if not (self._nnode is None):
            return self._nnode
        self._nnode = len(self.node_data)
        return self._nnode

           
    def write_gsf_file(self, filename="grid_spec.gsf"):
        """
        Write grid data to GSF file format.
        
        Args:
            filename (str): Output filename
        """
        try:
            filename = filename.replace(".txt", ".gsf")
            with open(filename, "w") as f:
                # Line 1: Header
                f.write("UNSTRUCTURED GWF\n")

                # Line 2: Global grid properties
                iz = 1  # Elevations of node and mesh element vertices 
                        # are supplied
                ic = 1  # Cell specifications associated with each node 
                        # are supplied
                f.write(f"{self.nnode} {self.nlays} {iz} {ic}\n")

                # Line 3: Number of vertices
                nvertex = len(self.vertices)
                f.write(f"{nvertex}\n")

                # Next NVERTEX lines: Vertex coordinates from DataFrame
                for _, row in self.vertices.iterrows():
                    f.write(f"{row['x']} {row['y']} {row['z']}\n")

                # Next NNODE lines: Node and element data from DataFrame
                for _, row in self.node_data.iterrows():
                    inode = row['inode']
                    x = row['x']
                    y = row['y']
                    z = row['z']
                    lay = row['lay']
                    vertex_indices = row['vertex_indices']
                    m = len(vertex_indices)

                    # Create a space-separated string of vertex indices
                    vertex_indices_str = " ".join(map(str, vertex_indices))

                    f.write(f"{inode} {x} {y} {z} {lay} {m} "
                           f"{vertex_indices_str}\n")

            print(f"Successfully wrote grid data to {filename}")

        except IOError as e:
            print(f"Error writing to file: {e}")
    
    @classmethod
    def read_gsf_file(cls, filename):
        """
        Read grid data from GSF file format and return a new Grid instance.
        
        Args:
            filename (str): Input filename
            
        Returns:
            Grid: New Grid instance with populated data
        """
        try:
            with open(filename, "r") as f:
                lines = f.readlines()
            
            # Remove newline characters and strip whitespace
            lines = [line.strip() for line in lines if line.strip()]
            
            # Line 1: Header
            if lines[0] != "UNSTRUCTURED GWF":
                raise ValueError("Invalid GSF file format: missing header")
            
            # Line 2: Global grid properties
            global_props = lines[1].split()
            if len(global_props) != 4:
                raise ValueError("Invalid GSF file format: incorrect global "
                                "properties line")
            
            nnode = int(global_props[0])
            nlay = int(global_props[1])
            iz = int(global_props[2])  # Elevations flag
            ic = int(global_props[3])  # Cell specifications flag
            
            # Line 3: Number of vertices
            nvertex = int(lines[2])
            
            # Read vertices (lines 4 to 4+nvertex-1)
            vertices_data = {'x': [], 'y': [], 'z': []}
            for i in range(3, 3 + nvertex):
                vertex_line = lines[i].split()
                if len(vertex_line) != 3:
                    raise ValueError(f"Invalid vertex line {i+1}: {lines[i]}")
                
                vertices_data['x'].append(float(vertex_line[0]))
                vertices_data['y'].append(float(vertex_line[1]))
                vertices_data['z'].append(float(vertex_line[2]))
            
            vertices_df = pd.DataFrame(vertices_data)
            
            # Read node data (lines 4+nvertex to end)
            node_data = {'inode': [], 'x': [], 'y': [], 'z': [], 'lay': [], 
                        'vertex_indices': []}
            
            for i in range(3 + nvertex, len(lines)):
                node_line = lines[i].split()
                if len(node_line) < 7:
                    raise ValueError(f"Invalid node line {i+1}: {lines[i]}")
                
                inode = int(node_line[0])
                x = float(node_line[1])
                y = float(node_line[2])
                z = float(node_line[3])
                lay = int(node_line[4])
                m = int(node_line[5])  # Number of vertices for this node
                
                # Read vertex indices
                if len(node_line) != 6 + m:
                    raise ValueError(f"Invalid node line {i+1}: expected "
                                    f"{6 + m} values, got {len(node_line)}")
                
                vertex_indices = [int(node_line[6 + j]) for j in range(m)]
                
                node_data['inode'].append(inode)
                node_data['x'].append(x)
                node_data['y'].append(y)
                node_data['z'].append(z)
                node_data['lay'].append(lay)
                node_data['vertex_indices'].append(vertex_indices)
            
            node_data_df = pd.DataFrame(node_data)
            
            # Verify we read the expected number of nodes
            if len(node_data_df) != nnode:
                raise ValueError(f"Mismatch: expected {nnode} nodes, "
                                f"read {len(node_data_df)}")
            
            # Create new Grid instance
            grid = cls()
            # grid.nnode = nnode
            grid.nlay = nlay
            grid.vertices = vertices_df
            grid.node_data = node_data_df
            
            print(f"Successfully read grid data from {filename}")
            print(f"Grid contains {grid.nnode} nodes, {grid.nlay} layers, "
                  f"{nvertex} vertices")
            
            return grid
            
        except IOError as e:
            print(f"Error reading file: {e}")
            raise
        except (ValueError, IndexError) as e:
            print(f"Error parsing GSF file: {e}")
            raise
        
    def _create_geometry(self):
        """
        Create a dataframe with columns node_id, geometry, top, botom
        """
        # loop through the node_data dataframe and create a shapely polygon
        # geometry for each node
        self._geometry = pd.DataFrame(columns=['inode', 'geometry', 'top', 
                                              'botm'])
        for _, row in self.node_data.iterrows():          
            node_id = row['inode']
            vertex_indices = row['vertex_indices']
            geometry = self.vertices.iloc[vertex_indices]
            shapely_ploygon = Polygon(geometry)
            self._geometry.loc[node_id, 'geometry'] = shapely_ploygon
            self._geometry.loc[node_id, 'top'] = geometry['z'].max()
            self._geometry.loc[node_id, 'botm'] = geometry['z'].min()

            # add centroid to geometry
            centroid = shapely_ploygon.centroid
            self._geometry.loc[node_id, 'xc'] = centroid.x
            self._geometry.loc[node_id, 'yc'] = centroid.y
            self._geometry.loc[node_id, 'zc'] = ((geometry['z'].max() + 
                                                 geometry['z'].min()) / 2)
            self._geometry.loc[node_id, 'inode'] = node_id

        self._geometry['layer'] = None
        for _, row in self._geometry.iterrows():
            if row['layer'] is not None:
                continue
            x = row['xc']
            y = row['yc']
            top = row['top']

            # get all nodes sharing (x,y)
            nodes_sharing_xy = self._geometry[(self._geometry['xc'] == x) & 
                                             (self._geometry['yc'] == y)]
            unique_tops = nodes_sharing_xy['top'].unique().tolist()
            unique_tops.sort(reverse=True)
            for i, top in enumerate(unique_tops):
                mask = nodes_sharing_xy['top'] == top
                self._geometry.loc[nodes_sharing_xy[mask].index, 'layer'] = i + 1

    
    @property
    def geometry(self):
        if self._geometry is None:
            self._create_geometry()
        return self._geometry
    

    def from_shapefile(self, filename: str, top, thickness) -> None:
        """
        # todo: test
        Create a grid from a shapefile. The shapefile is 2D with no z values.
        The top parameter specifies top elevation values, and thickness 
        specifies layer thicknesses.
        
        Args:
            filename: Path to the shapefile
            top: Column name (str), constant value (float/int), or None 
                 (default 1.0)
            thickness: List of column names (str), constant values 
                      (float/int), or None (default 1.0)
        """
        from typing import Union, List
        
        # Read the shapefile
        gdf = gpd.read_file(filename)
        
        # Handle top parameter
        if isinstance(top, str):
            # Column name in shapefile
            if top not in gdf.columns:
                raise ValueError(f"Column '{top}' not found in shapefile")
            top_values = gdf[top].values
        elif isinstance(top, (int, float)):
            # Constant value for all nodes
            top_values = [float(top)] * len(gdf)
        elif top is None:
            # Default constant value
            top_values = [1.0] * len(gdf)
        else:
            raise ValueError("top must be a string (column name), number "
                            "(constant), or None")
        
        # Handle thickness parameter
        if isinstance(thickness, list):
            # List of values
            thickness_values = []
            for i, thick in enumerate(thickness):
                if isinstance(thick, str):
                    # Column name in shapefile
                    if thick not in gdf.columns:
                        raise ValueError(f"Column '{thick}' not found in shapefile")
                    thickness_values.append(gdf[thick].values)
                elif isinstance(thick, (int, float)):
                    # Constant value for all nodes
                    thickness_values.append([float(thick)] * len(gdf))
                elif thick is None:
                    # Default constant value
                    thickness_values.append([1.0] * len(gdf))
                else:
                    raise ValueError(f"thickness[{i}] must be a string, "
                                    f"number, or None")       
        # Determine number of layers
        self.nlays = len(thickness_values)
        self.nnode = len(gdf) * self.nlays
        
        # Create vertices DataFrame
        vertices_data = []
        vertex_id = 0
        
        # Create vertices for each node and layer
        for node_idx in range(self.nnode):
            geom = gdf.iloc[node_idx].geometry
            top_val = top_values[node_idx]
            
            # Calculate bottom elevation for each layer
            current_top = top_val
            for layer_idx in range(self.nlays):
                thick_val = thickness_values[layer_idx][node_idx]
                bottom_val = current_top - thick_val
                
                # Create vertices for this layer (assuming rectangular cells)
                if hasattr(geom, 'exterior'):
                    # Polygon geometry
                    coords = list(geom.exterior.coords[:-1])  # Remove duplicate 
                                                              # last point
                else:
                    # Point geometry - create a small square around the point
                    x, y = geom.x, geom.y
                    size = 0.5  # Default cell size
                    coords = [
                        (x - size/2, y - size/2),
                        (x + size/2, y - size/2),
                        (x + size/2, y + size/2),
                        (x - size/2, y + size/2)
                    ]
                
                # Add vertices for top and bottom of this layer
                for x, y in coords:
                    # Top vertices
                    vertices_data.append({
                        'x': x,
                        'y': y,
                        'z': current_top
                    })
                    vertex_id += 1
                    
                    # Bottom vertices
                    vertices_data.append({
                        'x': x,
                        'y': y,
                        'z': bottom_val
                    })
                    vertex_id += 1
                
                current_top = bottom_val  # Next layer starts at bottom of 
                                          # current layer
        
        self.vertices = pd.DataFrame(vertices_data)
        
        # Create node_data DataFrame
        node_data = []
        vertex_id = 0
        
        for node_idx in range(self.nnode):
            geom = gdf.iloc[node_idx].geometry
            top_val = top_values[node_idx]
            
            # Calculate center coordinates
            if hasattr(geom, 'centroid'):
                center_x = geom.centroid.x
                center_y = geom.centroid.y
            else:
                center_x = geom.x
                center_y = geom.y
            
            current_top = top_val
            for layer_idx in range(self.nlays):
                thick_val = thickness_values[layer_idx][node_idx]
                bottom_val = current_top - thick_val
                center_z = (current_top + bottom_val) / 2
                
                # Calculate vertex indices for this cell
                if hasattr(geom, 'exterior'):
                    num_vertices = len(list(geom.exterior.coords)) - 1
                else:
                    num_vertices = 4  # Square around point
                
                # Each layer has 2 * num_vertices vertices (top + bottom)
                vertex_indices = list(range(vertex_id, vertex_id + 2 * num_vertices))
                vertex_id += 2 * num_vertices
                
                node_data.append({
                    'inode': node_idx * self.nlays + layer_idx,
                    'x': center_x,
                    'y': center_y,
                    'z': center_z,
                    'lay': layer_idx + 1,
                    'vertex_indices': vertex_indices
                })
                
                current_top = bottom_val
        
        self.node_data = pd.DataFrame(node_data)
        
        # Create geometry DataFrame
        geometry_data = []
        for node_idx in range(self.nnode):
            geom = gdf.iloc[node_idx].geometry
            top_val = top_values[node_idx]
            
            # Calculate total thickness
            total_thickness = sum(thickness_values[layer_idx][node_idx] 
                                for layer_idx in range(self.nlays))
            
            geometry_data.append({
                'node_id': node_idx,
                'geometry': geom,
                'top': top_val,
                'thickness': total_thickness
            })
        
        self._geometry = pd.DataFrame(geometry_data)

    def to_shapefile(self, filename: str) -> None:
        # todo: test
        # use data geometry to write shapefile using geopandas, also adds 
        # a column for node_id, top, thickness
        gdf = gpd.GeoDataFrame(self.geometry, geometry='geometry')
        gdf['node_id'] = self.geometry['node_id']
        gdf['top'] = self.geometry['top']
        gdf['thickness'] = self.geometry['thickness']
        gdf.to_file(filename)
    
    @classmethod
    def from_3d_polygon_list(cls, polygon_list: List[Polygon]) -> 'Grid':
        """
        Create a grid object from a list of polygons.        
        Args:
            polygon_list: List of Shapely Polygon objects            
        Returns:
            Grid: New grid instance with vertices and node data
        """
        # Initialize data structures
        vertices_data = []  # List to store unique vertices
        vertex_to_id = {}  # Map (x, y, z) -> vertex_id
        next_vertex_id = 0
        node_data = []  # List to store node information
        
        for ipolygon, polygon in enumerate(polygon_list):             
            current_vertices = list(polygon.exterior.coords[:-1])  # Remove last 
                                                                   # duplicate point
            vdf = pd.DataFrame(current_vertices, columns=['x', 'y', 'z'])
            unique_z = vdf['z'].unique()
            if len(unique_z) != 2:
                print(f"Warning: cell {ipolygon} is not a prism")
            else:
                ntop = np.sum(vdf['z'] == unique_z[0])
                nbot = np.sum(vdf['z'] == unique_z[1])
                if not (ntop == nbot):    
                    print(f"Warning: cell {ipolygon} is not a prism")
            vertex_indices = []  # Store vertex IDs for this polygon
            
            # Process each vertex in the current polygon
            
            for vertex_coords in current_vertices:
                x, y, z = vertex_coords                
               
                # Check if vertex already exists
                vertex_key = (x, y, z)
                if vertex_key not in vertex_to_id:
                    # Add new vertex
                    vertices_data.append({
                        'x': x,
                        'y': y, 
                        'z': z
                    })
                    vertex_to_id[vertex_key] = next_vertex_id
                    vertex_indices.append(next_vertex_id)
                    next_vertex_id += 1
                else:
                    # Use existing vertex
                    vertex_indices.append(vertex_to_id[vertex_key])
         
           
            # Add node data for this polygon
            xc,yc, zc = np.array(polygon.exterior.coords)[:-1].mean(axis = 0)
            node_data.append({
                'inode': ipolygon,  # 0-based indexing
                'x': xc,
                'y': yc,
               'z': zc,  # Defa0.0t z-coordinate
                'lay': None,  # Default layer
                'vertex_indices': vertex_indices
            })

        # Create grid instance
        grid = cls()               
        grid.vertices = pd.DataFrame(vertices_data)
        grid.node_data = pd.DataFrame(node_data)
       
        
        return grid

    def extrude(self, thickness):

        pass      

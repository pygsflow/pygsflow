from typing import List, Optional
import numpy as np
from shapely.geometry import Polygon
from .grid import Grid
import math
import pandas as pd

class StructuredGrid(Grid):
    def __init__(self, 
                 delr: np.ndarray, 
                 delc: np.ndarray, 
                 nrows: int, 
                 ncols: int,
                 nlays: int,
                 top: np.ndarray,  # 2D array (nrows, ncols)
                 botm: np.ndarray,  # 3D array (nlay, nrows, ncols)
                 origin_x: float = 0.0, 
                 origin_y: float = 0.0, 
                 rotation: float = 0.0,
                 filename: Optional[str] = None) -> None:
        """
        Initialize a structured grid.
        
        Args:
            delr: Array of grid cell sizes in x-direction
            delc: Array of grid cell sizes in y-direction
            nrows: Number of rows
            ncols: Number of columns
            top: 2D array of top elevations (nrows, ncols)
            botm: 3D array of bottom elevations (nlay, nrows, ncols)
            origin_x: X-coordinate of grid origin
            origin_y: Y-coordinate of grid origin
            rotation: Grid rotation angle in degrees
            filename: Optional filename for the grid
        """
        super().__init__()
        self.delr = delr
        self.delc = delc
        self.top = top
        self.botm = botm
        self.nrows = nrows
        self.ncols = ncols
        #self.nlays = botm.shape[0]
        self.origin_x = origin_x
        self.origin_y = origin_y
        self.rotation = rotation
        self.create()

    def create(self) -> List[Polygon]:
        """
        Creates a structured 3D grid of Shapely Polygons using cell size, 
        origin, and rotation. Creates 2D polygons representing the top 
        face of each 3D cell.
        
        Returns:
            A list of Shapely Polygon objects (top faces of 3D cells)
        """
        rows = self.nrows
        cols = self.ncols
        lays = self.botm.shape[0]
        total_cells = rows * cols * lays
        
        # Pre-compute rotation values (only once)
        if self.rotation == 0.0:
            # Fast path for no rotation
            self._create_no_rotation()
        
        # # todo: rotation is not tested
        # rotation_rad = math.radians(self.rotation)
        # cos_rot = math.cos(rotation_rad)
        # sin_rot = math.sin(rotation_rad)
        
        # # Pre-allocate list for better performance
        # grid: List[Polygon] = [None] * total_cells  # type: ignore
        
        # # Vectorized corner calculation
        # cell_indices = np.arange(total_cells)
        # row_indices = cell_indices // cols
        # col_indices = cell_indices % cols
        
        # # Calculate all cell corners at once (top face only for 3D grid)
        # x1 = col_indices * self.delr
        # y1 = row_indices * self.delc
        # x2 = (col_indices + 1) * self.delr
        # y2 = (row_indices + 1) * self.delc
        
        # # Create corner arrays for all cells (top face)
        # corners_x = np.array([x1, x2, x2, x1]).T  # (total_cells, 4)
        # corners_y = np.array([y1, y1, y2, y2]).T  # (total_cells, 4)
        
        # # Vectorized rotation and translation
        # x_rot = corners_x * cos_rot - corners_y * sin_rot
        # y_rot = corners_x * sin_rot + corners_y * cos_rot
        
        # x_final = x_rot + self.origin_x
        # y_final = y_rot + self.origin_y
        
        # # Create polygons efficiently (top faces of 3D cells)
        # for i in range(total_cells):
        #     corners = list(zip(x_final[i], y_final[i]))
        #     grid[i] = Polygon(corners)
        
        #return grid
    
    def _create_no_rotation(self) -> List[Polygon]:
        """
        Optimized creation for 3D grids with no rotation.
        Creates 2D polygons representing the top face of each 3D cell.
        
        Returns:
            A list of Shapely Polygon objects (top faces of 3D cells)
        """
        rows = self.nrows
        cols = self.ncols
        nlays = self.botm.shape[0]
        total_cells = rows * cols * nlays
        
        # Pre-allocate list
        grid = {} # type: ignore      
        vertix_dict = {} 
        node_data = []
        geometry = []
        v_id = 1
        for k in range(nlays):
            for i in range(rows):
                for j in range(cols):
                    cell_id = int(k * (rows * cols) + (i * cols) + j)
                    x1 = j * self.delr[j] + self.origin_x
                    y1 = i * self.delc[i] + self.origin_y
                    x2 = (j + 1) * self.delr[j] + self.origin_x
                    y2 = (i + 1) * self.delc[i] + self.origin_y
                    z1 = self.top[i, j] if k == 0 else self.botm[k-1, i, j]
                    z2 = self.botm[k, i, j]

                    xc = (x1 + x2) * 0.5
                    yc = (y1 + y2) * 0.5
                    zc = (z1 + z2) * 0.5

                    vertices = [
                        (x1, y1, z1),
                        (x2, y1, z1),
                        (x2, y2, z1),
                        (x1, y2, z1),
                        (x1, y1, z2),
                        (x2, y1, z2),
                        (x2, y2, z2),
                        (x1, y2, z2)
                    ]

                    node_vertex_list = []
                    for v in vertices:
                        if v not in vertix_dict.keys():
                            vertix_dict[v] = v_id
                            node_vertex_list.append(v_id)
                            v_id = v_id + 1
                        else:
                            node_vertex_list.append(vertix_dict[v])

                    grid[int(cell_id)] = Polygon(vertices)                    
                    node_data.append([cell_id+1, xc, yc, zc, k, node_vertex_list])
                    geometry.append([cell_id+1, grid[int(cell_id)], z1, z2])
        
        # Convert vertex dictionary to DataFrame
        vertices_data = []
        for (x, y, z), vertex_id in vertix_dict.items():
            vertices_data.append({'x': x, 'y': y, 'z': z, 
                                'vertex_id': vertex_id})
        self.vertices = pd.DataFrame(vertices_data)
        # Sort by vertex_id to maintain order
        self.vertices = self.vertices.sort_values('vertex_id').reset_index(
            drop=True)
        del(self.vertices['vertex_id'])
        
        self.node_data = pd.DataFrame(node_data, 
                                      columns=['inode', 'x', 'y', 'z', 'lay', 
                                              'vertex_indices'])

        self._geometry = pd.DataFrame(geometry, columns=['inode', 'geometry', 'top', 
                                              'botm'])
        
        
    # def _get_dis_vertices(self, dis):
    #     """
    #     Find unique vertices in a structured grid and compute cell information.
        
    #     Parameters:
    #     -----------
    #     dis : flopy discretization object
    #         Discretization package containing grid information
        
    #     Returns:
    #     --------
    #     vertices : numpy.ndarray
    #         Array of unique vertex coordinates (x, y, z)
    #     cell_info : dict
    #         Dictionary containing cell_id, vertex_lists, and centroids
    #     """
    #     delr = dis.delr.array
    #     delc = dis.delc.array
    #     nrow = dis.nrow
    #     ncol = dis.ncol
    #     nlay = dis.nlay
    #     top = dis.top.array
    #     bottom = dis.botm.array
        
    #     # Pre-calculate coordinate arrays
    #     x_coords = np.concatenate([[0], np.cumsum(delr)])
    #     y_coords = np.concatenate([[0], np.cumsum(delc)])
        
    #     # Initialize data structures
    #     vertices = []
    #     vertex_to_id = {}
    #     next_vertex_id = 1
        
    #     cell_ids = []
    #     vertex_lists = []
    #     centroids = []
    #     layers = []
    #     cols = []
    #     rows = []
        
    #     cell_id = 0
        
    #     # loop over layers (starting from layer 1 - top layer)
    #     for k in range(nlay):
    #         # loop over rows
    #         for i in range(nrow):
    #             # loop over columns
    #             for j in range(ncol):
    #                 # cell_id = cell_id + 1
    #                 cell_id = cell_id + 1
                    
    #                 # get the vertices (8 vertices). generate them in terms of 
    #                 # (i,j,k,"top") or (i,j,k,"botm")
    #                 # Calculate cell corner coordinates
    #                 x1, x2 = x_coords[j], x_coords[j+1]
    #                 y1, y2 = y_coords[i], y_coords[i+1]
                    
    #                 # Z coordinates for this cell
    #                 z_bottom = bottom[k, i, j]  # Bottom of current layer
    #                 if k == 0:
    #                     z_top = top[i, j]       # Top surface for layer 0
    #                 else:
    #                     z_top = bottom[k-1, i, j]  # Bottom of previous layer
                    
    #                 # 8 vertices of the cell (hexahedron)
    #                 vertex_coords = [
    #                     (x1, y1, z_top),     # 4: top-left-back
    #                     (x2, y1, z_top),     # 5: top-right-back
    #                     (x2, y2, z_top),     # 6: top-right-front
    #                     (x1, y2, z_top),      # 7: top-left-front
    #                     (x1, y1, z_bottom),  # 0: bottom-left-back
    #                     (x2, y1, z_bottom),  # 1: bottom-right-back
    #                     (x2, y2, z_bottom),  # 2: bottom-right-front
    #                     (x1, y2, z_bottom)  # 3: bottom-left-front
                        
    #                 ]
                    
    #                 # give each vertex a global unique id - start from 1.
    #                 vertex_ids = []
    #                 for coord in vertex_coords:
    #                     # before adding the vertices check if they already exist, 
    #                     # if they do, use the existing vertex id.
    #                     if coord not in vertex_to_id:
    #                         # add new vertices to the vertices list.
    #                         vertex_to_id[coord] = next_vertex_id
    #                         vertices.append(list(coord))
    #                         next_vertex_id += 1
    #                     vertex_ids.append(vertex_to_id[coord])
                    
    #                 # Calculate centroid
    #                 centroid = [
    #                     (x1 + x2) * 0.5,              # x
    #                     (y1 + y2) * 0.5,              # y
    #                     (z_top + z_bottom) * 0.5      # z
    #                 ]
                    
    #                 # add the cell_id to the cell_info dictionary.
    #                 cell_ids.append(cell_id)
    #                 # add the vertex_lists to the cell_info dictionary.
    #                 vertex_lists.append(vertex_ids)
    #                 # add the centroids to the cell_info dictionary.
    #                 centroids.append(centroid)
    #                 # add the layer to the cell_info dictionary.
    #                 layers.append(k)
    #                 cols.append(j)
    #                 rows.append(i)
        
    #     # Create cell_info dictionary
    #     cell_info = {
    #         'cell_ids': np.array(cell_ids),
    #         'vertex_lists': np.array(vertex_lists),
    #         'centroids': np.array(centroids),
    #         'layer': np.array(layers)
    #     }
        
    #     # Convert vertices to numpy array
    #     vertices = np.array(vertices)
        
    #     return vertices, cell_info     

    @classmethod
    def from_structured_dis(cls, dis):       
        # Create grid instance with required arguments
        # Handle delr and delc - convert scalars to constant arrays
        if np.isscalar(dis.delr):
            delr_array = np.full(dis.ncol, dis.delr)
        else:
            delr_array = dis.delr.array
            
        if np.isscalar(dis.delc):
            delc_array = np.full(dis.nrow, dis.delc)
        else:
            delc_array = dis.delc.array
        
        grid = cls(
            delr=delr_array,
            delc=delc_array,
            nrows=dis.nrow,
            ncols=dis.ncol,
            nlays=dis.nlay,
            top=getattr(dis.top, 'array', 0.0),
            botm=getattr(dis.botm, 'array', 0.0),
            origin_x=getattr(dis, 'origin_x', 0.0),
            origin_y=getattr(dis, 'origin_y', 0.0),
            rotation=getattr(dis, 'rotation', 0.0)
        )
        # grid.nlays = dis.nlay

        # grid.nnode = dis.nrow * dis.ncol * dis.nlay
        # grid.nlay = dis.nlay

        # # Get vertices and cell information
        # vertices, cell_info = grid._get_dis_vertices(dis)
       
        # grid.vertices = pd.DataFrame(vertices, columns=['x', 'y', 'z'])

        # # Create node data DataFrame
        # cell_nodes = cell_info['cell_ids']
        # node_data = []
        # for i, inode in enumerate(cell_nodes):
        #     x = cell_info['centroids'][i, 0]
        #     y = cell_info['centroids'][i, 1]
        #     z = cell_info['centroids'][i, 2]
        #     lay = cell_info['layer'][i]
        #     vertex_indices = cell_info['vertex_lists'][i]
        #     node_data.append([inode, x, y, z, lay, vertex_indices])
        
        # df_node_data = pd.DataFrame(node_data, 
        #                            columns=['inode', 'x', 'y', 'z', 'lay', 
        #                                    'vertex_indices'])
        # grid.node_data = df_node_data    

        return grid

    
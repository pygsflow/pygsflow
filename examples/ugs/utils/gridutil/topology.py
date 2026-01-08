# give a grid class as defined in grid_io.py, write a function that calculates the topology of the grid.

from typing import List, Tuple, Optional
from .grid import Grid
import pandas as pd
import numpy as np
from shapely.geometry import Polygon, Point
from shapely.ops import transform

# in the cell topology, always assume top face is identical to bottom face and parallel to the horizontal face.
class Topology():
    def __init__(self, grid: Grid) -> None:
        """
        Initialize GridTopology with a Grid instance.
        
        Args:
            grid: Grid instance with vertices and node_data DataFrames
        """
        self.grid = grid
        self._validate_grid()
        self._vertex_cache = {}  # Cache for vertex lookups
        self._cell_cache = {}    # Cache for cell data lookups
        self._face_cache = {}    # Cache for face-sharing results
        
        # Pre-compute frequently used data structures
        self._node_data_dict = self.grid.node_data.set_index('inode').to_dict('index')
        self._vertices_array = self.grid.vertices[['x', 'y', 'z']].values
    
    def _validate_grid(self) -> None:
        """Validate that the grid has the required attributes and data."""
        if not hasattr(self.grid, 'vertices') or not hasattr(self.grid, 'node_data'):
            raise ValueError("Grid must have 'vertices' and 'node_data' attributes")
        
        if not isinstance(self.grid.vertices, pd.DataFrame):
            raise ValueError("Grid.vertices must be a pandas DataFrame")
        
        if not isinstance(self.grid.node_data, pd.DataFrame):
            raise ValueError("Grid.node_data must be a pandas DataFrame")
        
        if self.grid.vertices.empty or self.grid.node_data.empty:
            raise ValueError("Grid vertices and node_data cannot be empty")
        
        # Validate vertices DataFrame columns
        expected_vertex_columns = ['x', 'y', 'z']
        if not all(col in self.grid.vertices.columns for col in expected_vertex_columns):
            missing_cols = [col for col in expected_vertex_columns if col not in self.grid.vertices.columns]
            raise ValueError(f"Grid.vertices missing required columns: {missing_cols}")
        
        # Validate node_data DataFrame columns
        expected_node_columns = ['inode', 'x', 'y', 'z', 'lay', 'vertex_indices']
        if not all(col in self.grid.node_data.columns for col in expected_node_columns):
            missing_cols = [col for col in expected_node_columns if col not in self.grid.node_data.columns]
            raise ValueError(f"Grid.node_data missing required columns: {missing_cols}")
        
        # Validate vertex_indices column contains lists
        if not all(isinstance(idx_list, list) for idx_list in self.grid.node_data['vertex_indices']):
            raise ValueError("Grid.node_data['vertex_indices'] must contain lists")
        
        # Validate that vertex indices are within valid range
        max_vertex_id = len(self.grid.vertices)
        for idx, row in self.grid.node_data.iterrows():
            for vertex_id in row['vertex_indices']:
                if not isinstance(vertex_id, int) or vertex_id < 1 or vertex_id > max_vertex_id:
                    raise ValueError(f"Invalid vertex ID {vertex_id} in cell {row['inode']}. "
                                   f"Vertex IDs must be integers between 1 and {max_vertex_id}")
    
    def _cells_sharing_vertex(self, vertex_id: int) -> pd.DataFrame:
        """
        Find all cells that share a given vertex.
        Optimized version with caching, early validation, and vectorized operations.
        
        Args:
            vertex_id: Vertex ID (1-based indexing)
            
        Returns:
            DataFrame of cells containing the vertex
        """
        # Check cache first
        if vertex_id in self._vertex_cache:
            return self._vertex_cache[vertex_id]
        
        # Early validation - check if vertex_id is in valid range
        if vertex_id < 1 or vertex_id > len(self.grid.vertices):
            empty_df = pd.DataFrame(columns=self.grid.node_data.columns)
            self._vertex_cache[vertex_id] = empty_df
            return empty_df
        
        # Use vectorized operations instead of apply with lambda
        # Convert vertex_indices to numpy arrays for faster operations
        vertex_arrays = self.grid.node_data['vertex_indices'].values
        
        # Create boolean mask using list comprehension (faster than apply)
        mask = [vertex_id in vertex_list for vertex_list in vertex_arrays]
        
        result = self.grid.node_data[mask]
        
        # Cache the result
        self._vertex_cache[vertex_id] = result
        
        return result
    
    def clear_vertex_cache(self) -> None:
        """
        Clear the vertex lookup cache.
        Useful when grid data changes.
        """
        self._vertex_cache.clear()
    
    def clear_all_caches(self) -> None:
        """
        Clear all caches (vertex, cell, and face caches).
        Useful when grid data changes.
        """
        self._vertex_cache.clear()
        self._cell_cache.clear()
        self._face_cache.clear()
    
    def get_cell_neighbors(self, cell_id: int) -> List[int]:
        """
        Find neighboring cells that share vertices with the given cell.
        Optimized version using pre-computed data structures.
        
        Args:
            cell_id: Cell ID (1-based indexing)
            
        Returns:
            List of neighboring cell IDs
        """
        # Check cache first
        if cell_id in self._cell_cache:
            return self._cell_cache[cell_id]
        
        # Get the cell's vertex indices using pre-computed dictionary
        if cell_id not in self._node_data_dict:
            raise ValueError(f"Cell {cell_id} not found")
        
        vertex_indices = self._node_data_dict[cell_id]['vertex_indices']
        
        # Find all cells that share any vertex with this cell
        neighbors = set()
        for vertex_id in vertex_indices:
            cells_with_vertex = self._cells_sharing_vertex(vertex_id)
            neighbor_ids = cells_with_vertex['inode'].tolist()
            neighbors.update(neighbor_ids)
        
        # Remove the original cell from neighbors
        neighbors.discard(cell_id)
        result = list(neighbors)
        
        # Cache the result
        self._cell_cache[cell_id] = result
        return result
    
    def get_vertex_coordinates(self, vertex_id: int) -> Tuple[float, float, float]:
        """
        Get coordinates of a specific vertex.
        Optimized version using pre-computed array.
        
        Args:
            vertex_id: Vertex ID (1-based indexing)
            
        Returns:
            (x, y, z) coordinates
        """
        if vertex_id < 1 or vertex_id > len(self.grid.vertices):
            raise ValueError(f"Vertex ID {vertex_id} out of range")
        
        # Use pre-computed array for faster access
        coords = self._vertices_array[vertex_id - 1]  # Convert to 0-based indexing
        return (coords[0], coords[1], coords[2])
    
    def get_cell_vertices(self, cell_id: int) -> pd.DataFrame:
        """
        Get all vertices of a specific cell.
        Optimized version using pre-computed data structures.
        
        Args:
            cell_id: Cell ID (1-based indexing)
            
        Returns:
            DataFrame with vertex coordinates
        """
        # Use pre-computed dictionary for faster access
        if cell_id not in self._node_data_dict:
            raise ValueError(f"Cell {cell_id} not found")
        
        vertex_indices = self._node_data_dict[cell_id]['vertex_indices']
        
        # Vectorized coordinate retrieval using pre-computed array
        vertex_coords = []
        for vid in vertex_indices:
            coords = self._vertices_array[vid - 1]  # Convert to 0-based indexing
            vertex_coords.append({'vertex_id': vid, 'x': coords[0], 'y': coords[1], 'z': coords[2]})
        
        return pd.DataFrame(vertex_coords)
    
    
    def calculate_cell_volume(self, cell_id: int) -> float:
        """
        Calculate volume of a cell using the simplified method:
        Volume = Top face area × Cell thickness
        
        Since top and bottom faces are identical and parallel, we can use
        the area of the top face multiplied by the cell thickness.
        
        Args:
            cell_id: Cell ID (1-based indexing)
            
        Returns:
            Cell volume
        """
        # Get all vertices for the cell
        vertices_df = self.get_cell_vertices(cell_id)
        
        # Calculate cell thickness (z-difference between top and bottom faces)
        thickness = self._calculate_cell_thickness(vertices_df)
        
        # Calculate top face area
        top_area = self._calculate_horizontal_face_area(cell_id, face_type='top')
        
        # Volume = top face area × thickness
        volume = top_area * thickness
        return volume
    
    def _calculate_horizontal_face_area(self, cell_id: int, face_type: str = 'top') -> float:
        """
        Calculate area of the horizontal face (top or bottom) of a cell.
        
        Args:
            cell_id: Cell ID (1-based indexing)
            face_type: 'top' or 'bottom' face
            
        Returns:
            Horizontal face area
        """
        # Get all vertices for the cell
        vertices_df = self.get_cell_vertices(cell_id)
        num_vertices = len(vertices_df)
        
        # Determine which vertices belong to the top or bottom face
        if face_type == 'top':
            # For top face, use vertices with highest z-coordinates
            z_coords = vertices_df['z'].values
            max_z = z_coords.max()
            face_vertices = vertices_df[vertices_df['z'] == max_z]
            face_z_value = max_z
        elif face_type == 'bottom':
            # For bottom face, use vertices with lowest z-coordinates
            z_coords = vertices_df['z'].values
            min_z = z_coords.min()
            face_vertices = vertices_df[vertices_df['z'] == min_z]
            face_z_value = min_z
        else:
            raise ValueError("face_type must be 'top' or 'bottom'")
        
        # Validate that all vertices in the face have identical z-coordinates
        face_z_coords = face_vertices['z'].values
        if not np.allclose(face_z_coords, face_z_value, rtol=1e-10, atol=1e-10):
            raise ValueError(f"Cell {cell_id} {face_type} face vertices do not have identical z-coordinates. "
                           f"Expected all z = {face_z_value}, but got: {face_z_coords}")
        
        # Calculate area using shoelace formula
        if len(face_vertices) < 3:
            raise ValueError(f"Face must have at least 3 vertices, got {len(face_vertices)}")
        
        # Use Shapely for optimized area calculation
        ordered_vertices = self._order_vertices_for_polygon(face_vertices)
        
        # Create Shapely polygon from ordered vertices
        coords = list(zip(ordered_vertices['x'], ordered_vertices['y']))
        polygon = Polygon(coords)
        
        # Get area using Shapely (much faster and more accurate)
        return polygon.area    
    
    def _order_vertices_for_polygon(self, face_vertices: pd.DataFrame) -> pd.DataFrame:
        """
        Order vertices of a horizontal face in counterclockwise order for proper polygon area calculation.
        
        This method ensures that vertices are ordered correctly for the shoelace formula
        by sorting them in counterclockwise order around the polygon's centroid.
        
        Args:
            face_vertices: DataFrame with vertices of the face (x, y, z columns)
            
        Returns:
            DataFrame with vertices ordered counterclockwise
        """
        if len(face_vertices) < 3:
            return face_vertices
        
        # Calculate centroid of the face
        centroid_x = face_vertices['x'].mean()
        centroid_y = face_vertices['y'].mean()
        
        # Calculate angles from centroid to each vertex
        angles = []
        for _, vertex in face_vertices.iterrows():
            # Calculate angle from centroid to vertex
            dx = vertex['x'] - centroid_x
            dy = vertex['y'] - centroid_y
            angle = np.arctan2(dy, dx)
            angles.append(angle)
        
        # Sort vertices by angle (counterclockwise order)
        sorted_indices = np.argsort(angles)
        ordered_vertices = face_vertices.iloc[sorted_indices].reset_index(drop=True)
        
        return ordered_vertices

    def _calculate_cell_thickness(self, vertices_df: pd.DataFrame) -> float:
        """
        Calculate the thickness of a cell (z-difference between top and bottom faces).
        
        Args:
            vertices_df: DataFrame with all cell vertices
            
        Returns:
            Cell thickness
        """
        z_coords = vertices_df['z'].values
        
        # Validate that we have exactly two distinct z-levels (top and bottom)
        unique_z = np.unique(z_coords)
        if len(unique_z) != 2:
            raise ValueError(f"Cell must have exactly 2 distinct z-levels (top and bottom), "
                           f"but found {len(unique_z)} levels: {unique_z}")
        
        # Validate that vertices are properly distributed between top and bottom
        min_z, max_z = unique_z[0], unique_z[1]        
        
        # Calculate thickness
        thickness = max_z - min_z
        if thickness <= 0:
            raise ValueError(f"Cell thickness must be positive, got: {thickness}")
        
        return thickness

    # def _cells_sharing_face(self, cell_id: int) -> pd.DataFrame:
    #     """
    #     Find all cells that share a face with the given cell.
    #     Optimized version with caching and pre-computed data structures.
        
    #     Args:
    #         cell_id: Cell ID (1-based indexing)
            
    #     Returns:
    #         DataFrame of cells that share a face with the given cell
    #     """
    #     # Check cache first
    #     if cell_id in self._face_cache:
    #         return self._face_cache[cell_id]
        
    #     # Get all vertex-sharing neighbors using existing function
    #     all_neighbors = self.get_cell_neighbors(cell_id)
        
    #     # Get the cell's vertex indices using pre-computed dictionary
    #     if cell_id not in self._node_data_dict:
    #         raise ValueError(f"Cell {cell_id} not found")
        
    #     vertex_indices = set(self._node_data_dict[cell_id]['vertex_indices'])
        
    #     # Filter to only cells that share enough vertices to form a face
    #     face_neighbors = []
    #     for neighbor_id in all_neighbors:
    #         neighbor_vertices = set(self._node_data_dict[neighbor_id]['vertex_indices'])
            
    #         # Calculate shared vertices
    #         shared_vertices = vertex_indices.intersection(neighbor_vertices)
            
    #         # For face-sharing, cells must share at least 3 vertices (minimum for a face)
    #         if len(shared_vertices) >= 3:
    #             # Validate that shared vertices form a valid polygon face
    #             if self._validate_shared_vertices_form_face(cell_id, neighbor_id, shared_vertices):
    #                 face_neighbors.append(neighbor_id)
        
    #     # Return data as follows: cell_id, neighbor_id, area, distance
    #     if not face_neighbors:
    #         data = pd.DataFrame(columns=['cell_id', 'neighbor_id', 'area', 'distance'])
    #     else:
    #         # Build data efficiently using list comprehension
    #         data_rows = []
    #         for neighbor_id in face_neighbors:
    #             # Use optimized Shapely methods
    #             area = self._calculate_face_area_shapely(cell_id, neighbor_id)
    #             distance = self._calculate_distance_to_face(cell_id, neighbor_id)
    #             data_rows.append({
    #                 'cell_id': cell_id, 
    #                 'neighbor_id': neighbor_id, 
    #                 'area': area, 
    #                 'distance': distance
    #             })
    #         data = pd.DataFrame(data_rows)
        
    #     # Cache the result
    #     self._face_cache[cell_id] = data
    #     return data
        

    
    def _validate_shared_vertices_form_face(self, cell_id: int, neighbor_id: int, shared_vertices: set) -> bool:
        """
        Validate that shared vertices between two cells form a valid polygon face.
        Optimized version with early returns and vectorized operations.
        
        Args:
            cell_id: First cell ID
            neighbor_id: Second cell ID  
            shared_vertices: Set of shared vertex IDs
            
        Returns:
            True if shared vertices form a valid face, False otherwise
        """
        if len(shared_vertices) < 3:
            return False
        
        # Get coordinates of shared vertices using pre-computed array
        vertex_ids = list(shared_vertices)
        coords = self._vertices_array[[vid - 1 for vid in vertex_ids]]
        
        # Early return: Check if vertices are degenerate (all same point)
        if self._are_vertices_degenerate_fast(coords):
            return False
        
        # Early return: Check if vertices are collinear
        if self._are_vertices_collinear_fast(coords):
            return False
        
        # Check if vertices form a polygon with non-zero area
        return self._has_nonzero_area_fast(coords)
    
    def _are_vertices_collinear_fast(self, coords: np.ndarray) -> bool:
        """
        Fast check if vertices are collinear (all on the same line).
        Optimized version using vectorized operations.
        
        Args:
            coords: Array of vertex coordinates (n x 3)
            
        Returns:
            True if vertices are collinear, False otherwise
        """
        if len(coords) < 3:
            return True  # Less than 3 points are always collinear
        
        # For 3 points, check if they are collinear
        if len(coords) == 3:
            v1 = coords[1] - coords[0]
            v2 = coords[2] - coords[0]
            cross_product = np.cross(v1, v2)
            return np.allclose(cross_product, 0, atol=1e-10)
        
        # For 4+ points, check if all points lie on the line defined by first 2 points
        p1, p2 = coords[0], coords[1]
        direction = p2 - p1
        
        # If direction vector is zero, all points are the same
        if np.allclose(direction, 0, atol=1e-10):
            return True
        
        # Check if all other points lie on the line
        for i in range(2, len(coords)):
            point = coords[i]
            # Vector from p1 to current point
            v = point - p1
            # Cross product should be zero for collinear points
            cross_product = np.cross(direction, v)
            if not np.allclose(cross_product, 0, atol=1e-10):
                return False
        
        return True
    
    def _has_nonzero_area_fast(self, coords: np.ndarray) -> bool:
        """
        Fast check if vertices form a polygon with non-zero area.
        Optimized version using Shapely.
        
        Args:
            coords: Array of vertex coordinates (n x 3)
            
        Returns:
            True if polygon has non-zero area, False otherwise
        """
        if len(coords) < 3:
            return False
        
        # Use the plane with the largest projection area
        # Check xy, xz, and yz projections and use the one with largest area
        areas = []
        
        # XY projection using Shapely
        try:
            xy_coords = [(x, y) for x, y in zip(coords[:, 0], coords[:, 1])]
            xy_polygon = Polygon(xy_coords)
            areas.append(xy_polygon.area)
        except:
            areas.append(0.0)
        
        # XZ projection using Shapely
        try:
            xz_coords = [(x, z) for x, z in zip(coords[:, 0], coords[:, 2])]
            xz_polygon = Polygon(xz_coords)
            areas.append(xz_polygon.area)
        except:
            areas.append(0.0)
        
        # YZ projection using Shapely
        try:
            yz_coords = [(y, z) for y, z in zip(coords[:, 1], coords[:, 2])]
            yz_polygon = Polygon(yz_coords)
            areas.append(yz_polygon.area)
        except:
            areas.append(0.0)
        
        # Use the projection with the largest area
        max_area = max(areas)
        return max_area > 1e-10  # Non-zero area threshold
    
    def _are_vertices_degenerate_fast(self, coords: np.ndarray) -> bool:
        """
        Fast check if vertices are degenerate (all the same point or very close).
        Optimized version using vectorized operations.
        
        Args:
            coords: Array of vertex coordinates (n x 3)
            
        Returns:
            True if vertices are degenerate, False otherwise
        """
        if len(coords) < 2:
            return True
        
        # Vectorized check: compute distances from first point to all others
        first_point = coords[0]
        distances = np.linalg.norm(coords - first_point, axis=1)
        
        # Check if all distances are within tolerance
        return np.all(distances < 1e-10)
    
    def _calculate_distance_to_face(self, cell_id: int, neighbor_id: int) -> float:
        """
        Calculate the distance between the centers of two face-sharing cells.
        Optimized version using pre-computed data structures and Shapely.
        
        Args:
            cell_id: First cell ID
            neighbor_id: Second cell ID
            
        Returns:
            Distance between cell centers
        """
        # Get cell centers using pre-computed dictionary
        if cell_id not in self._node_data_dict or neighbor_id not in self._node_data_dict:
            raise ValueError(f"Cell {cell_id} or {neighbor_id} not found")
        
        cell_data = self._node_data_dict[cell_id]
        neighbor_data = self._node_data_dict[neighbor_id]
        
        # Get cell center coordinates
        cell_center = Point(cell_data['x'], cell_data['y'], cell_data['z'])
        neighbor_center = Point(neighbor_data['x'], neighbor_data['y'], neighbor_data['z'])
        
        # Calculate distance using Shapely
        return cell_center.distance(neighbor_center)
    
    def _calculate_face_area_shapely(self, cell_id: int, neighbor_id: int) -> float:
        """
        Calculate the area of the shared face between two cells using Shapely.
        Optimized version using pre-computed data structures.
        
        Args:
            cell_id: First cell ID
            neighbor_id: Second cell ID
            
        Returns:
            Area of the shared face
        """
        # Get shared vertices between the two cells using pre-computed dictionary
        if cell_id not in self._node_data_dict or neighbor_id not in self._node_data_dict:
            raise ValueError(f"Cell {cell_id} or {neighbor_id} not found")
        
        cell_vertices = set(self._node_data_dict[cell_id]['vertex_indices'])
        neighbor_vertices = set(self._node_data_dict[neighbor_id]['vertex_indices'])
        shared_vertices = cell_vertices.intersection(neighbor_vertices)
        
        if len(shared_vertices) < 3:
            return 0.0
        
        # Get coordinates of shared vertices using pre-computed array
        vertex_ids = list(shared_vertices)
        coords = self._vertices_array[[vid - 1 for vid in vertex_ids]]
        
        # Step 1: Fast coplanarity check
        if not self._validate_coplanar_polygon_fast(coords):
            return 0.0
        
        # Step 2: Order vertices around centroid for accurate area calculation
        ordered_coords = self._order_vertices_around_centroid(coords)
        
        # Step 3: Calculate area using fast cross product method
        return self._calculate_coplanar_polygon_area_fast(ordered_coords)
    
    def _validate_coplanar_polygon_fast(self, coords: np.ndarray) -> bool:
        """
        Fast coplanarity check using vector mathematics.
        
        Args:
            coords: Array of vertex coordinates (n_vertices, 3)
            
        Returns:
            True if vertices are coplanar, False otherwise
        """
        if len(coords) < 3:
            return False
        
        # Fast path for 3 vertices (always coplanar)
        if len(coords) == 3:
            return True
        
        # Use first 3 points to define the plane
        p1, p2, p3 = coords[0], coords[1], coords[2]
        
        # Calculate plane normal vector
        v1 = p2 - p1
        v2 = p3 - p1
        normal = np.cross(v1, v2)
        
        # Check if normal vector is zero (collinear points)
        normal_norm = np.linalg.norm(normal)
        if normal_norm < 1e-10:
            return False
        
        # Normalize normal vector for efficiency
        normal = normal / normal_norm
        
        # Check that all other points lie on the plane (vectorized)
        if len(coords) > 3:
            other_points = coords[3:]
            vectors = other_points - p1
            # Dot product with normal should be ~0 for coplanar points
            distances = np.abs(np.dot(vectors, normal))
            return np.all(distances < 1e-10)
        
        return True
    
    def _order_vertices_around_centroid(self, coords: np.ndarray) -> np.ndarray:
        """
        Order vertices counterclockwise around the centroid using proper 2D plane projection.
        
        Args:
            coords: Array of vertex coordinates (n_vertices, 3)
            
        Returns:
            Ordered array of vertex coordinates
        """
        if len(coords) < 3:
            return coords
        
        # Step 1: Find the centroid
        centroid = np.mean(coords, axis=0)
        
        # Step 2: Project to a 2D plane using the polygon's plane basis vectors
        # Use the same vectors from coplanarity check to define the plane
        p1, p2, p3 = coords[0], coords[1], coords[2]
        v1 = p2 - p1  # First basis vector
        v2 = p3 - p1  # Second basis vector
        
        # Normalize the basis vectors
        v1_norm = v1 / np.linalg.norm(v1)
        v2_norm = v2 / np.linalg.norm(v2)
        
        # Make v2 orthogonal to v1 using Gram-Schmidt process
        v2_ortho = v2_norm - np.dot(v2_norm, v1_norm) * v1_norm
        v2_ortho = v2_ortho / np.linalg.norm(v2_ortho)
        
        # Step 3: Project all points to the 2D plane
        # Translate points relative to centroid
        translated_coords = coords - centroid
        
        # Project onto the 2D plane using the basis vectors
        # Each point becomes (dot(v1_norm, point), dot(v2_ortho, point))
        proj_coords = np.column_stack([
            np.dot(translated_coords, v1_norm),
            np.dot(translated_coords, v2_ortho)
        ])
        
        # Step 4: Perform angular sort
        # Calculate angles from centroid (which is now at origin) to each vertex
        angles = np.arctan2(proj_coords[:, 1], proj_coords[:, 0])
        
        # Sort by angle (counterclockwise from positive x-axis)
        sorted_indices = np.argsort(angles)
        return coords[sorted_indices]
    
    def _calculate_coplanar_polygon_area_fast(self, coords: np.ndarray) -> float:
        """
        Fast area calculation using ordered vertices and cross product method.
        
        Args:
            coords: Array of vertex coordinates (n_vertices, 3) - must be coplanar and ordered
            
        Returns:
            True area of the polygon face in 3D space
        """
        if len(coords) < 3:
            return 0.0
        
        # Use the optimized cross product method for 3D polygon area
        # This is the standard mathematical approach for 3D polygon area calculation
        
        total = np.array([0.0, 0.0, 0.0])
        
        # Sum cross products of consecutive vertices (vectorized where possible)
        for i in range(len(coords)):
            vi1 = coords[i]
            vi2 = coords[(i + 1) % len(coords)]  # Wrap around to first vertex
            total += np.cross(vi1, vi2)
        
        # Area is half the magnitude of the resulting vector
        return np.linalg.norm(total) / 2.0
    

    def grid_topology(self, grid: Grid) -> None:
        """
        Calculate the topology of the grid.
        """
        self._validate_grid()

        # creat a dataframe to store the topology of the grid. the first column is cell_id
        # and the second column is the neighbor_id, the thrid column is volume of the cell, the
        # fourth column is the area of the face, the fifth column is the distance between the 
        # centers of the two cells. the function must be optimal cell that topolgy of cell_id = 1 and cell_id = 2
        # is the same as the topology of cell_id = 2 and cell_id = 1.
        topology_data = []  # List to collect all topology data
        processed_pairs = set()  # Track processed cell pairs to avoid duplicates
        
        for cell_id in range(1, self.grid.nnode + 1):
            neighbors = self.get_cell_neighbors(cell_id)
            # Calculate volume once per cell (not per neighbor)
            volume = self.calculate_cell_volume(cell_id)
            
            for neighbor_id in neighbors:
                # Create a unique pair identifier (smaller_id, larger_id)
                pair = (min(cell_id, neighbor_id), max(cell_id, neighbor_id))
                
                # Only process if we haven't seen this pair before
                if pair not in processed_pairs:
                    processed_pairs.add(pair)
                    
                    # Calculate area and distance once per unique pair
                    area = self._calculate_face_area_shapely(cell_id, neighbor_id)
                    distance = self._calculate_distance_to_face(cell_id, neighbor_id)
                    
                    # Get neighbor's volume for the pair
                    neighbor_volume = self.calculate_cell_volume(neighbor_id)
                    
                    # Add both directions of the pair to the data list
                    if area > 0:
                        topology_data.extend([
                            {
                                'cell_id': cell_id, 
                                'neighbor_id': neighbor_id, 
                                'volume': volume, 
                                'area': area, 
                                'distance': distance
                            },
                            {
                                'cell_id': neighbor_id, 
                                'neighbor_id': cell_id, 
                                'volume': neighbor_volume, 
                                'area': area, 
                                'distance': distance
                            }
                        ])
        # Create DataFrame from collected data
        topology_df = pd.DataFrame(topology_data)
        return topology_df

       

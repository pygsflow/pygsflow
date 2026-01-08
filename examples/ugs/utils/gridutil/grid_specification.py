import pandas as pd


"""

This documentation is extracted from PEST Groundwater Data Utilities
Part A, page 21

File Specification:
# ---------------------------------------------------------------------
Line 1: "UNSTRUCTURED GWF"
Note: If "GWF" is omitted, then it is assumed.

Line 2: nnode nlay iz ic
where:
	- nnode is the number of nodes in the grid
	nlay is the number of layers in the model
	- iz is 1 if elevations of node and mesh element vertices are
	supplied; 0 otherwise
	- ic is 1 if the cell specifications associated with each node are
	supplied; 0 otherwise

	Note. All utilities documented herein, as well as PLPROC, require
	that both iz and ic be 1. Options associated with values other
	than 1 are not presented below. If these are omitted from the grid
	specification file they are assumed to be 1.

A list of vertex coordinates is now provided. These vertices will
	be cited by index later in this file where mesh geometric details
	are provided. The indices of vertices are defined implicitly by
	their positions in the following list. Numbering is assumed to
	begin at 1.
	
Line 3: nvertex 
where:
	nvertex is the number of element vertex definitions to follow
	
NVERTEX next lines:
x, y, z
where:
	x, y and z are the coordinates of each vertex

NNODE next lines: inode x y z lay m (ivertex(i),i=1,m)
where:
	- inode is a node number (these must be supplied in increasing
	order starting from 1)
	- x, y and z are node coordinates
	- lay is the layer number of the node
	- m is the number of vertices defining the three-dimensional
	element associated with the node
	- ivertex(i) are vertex indices defined with reference to the vertex
	list provided above
"""
def read_gsf_file(filename):
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
            if lines[0][0]== "#":
                del(lines[0])

            
            # Line 1: Header
            if not("UNSTRUCTURED"  in lines[0]):
                raise ValueError(
                    "Invalid GSF file format: missing header")
            
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
                    raise ValueError(
                        f"Invalid vertex line {i+1}: {lines[i]}")
                
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
            # class Grid:
            #     pass
            # grid = Grid()
            # # grid.nnode = nnode
            # grid.nlay = nlay
            # grid.vertices = vertices_df
            # grid.node_data = node_data_df
            
            print(
                f"Successfully read grid data from {filename}")
            print(
                f"Grid contains {nnode} nodes, {nlay} layers, "
                f"{nvertex} vertices")
            
            return vertices_df, node_data_df
            
        except IOError as e:
            print(f"Error reading file: {e}")
            raise
        except (ValueError, IndexError) as e:
            print(f"Error parsing GSF file: {e}")
            raise
    
def write_gsf_file(vertices, node_data, filename="grid_spec.gsf"):
    """
    Write grid data to GSF file format.
    
    Args:
        filename (str): Output filename
    """
    nnode = len(node_data)
    nlays = node_data['lay'].max()
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
            f.write(f"{nnode} {nlays} {iz} {ic}\n")

            # Line 3: Number of vertices
            nvertex = len(vertices)
            f.write(f"{nvertex}\n")

            # Next NVERTEX lines: Vertex coordinates from
            # DataFrame
            for _, row in vertices.iterrows():
                f.write(f"{row['x']} {row['y']} {row['z']}\n")

            # Next NNODE lines: Node and element data from
            # DataFrame
            for _, row in node_data.iterrows():
                inode = row['inode']
                x = row['x']
                y = row['y']
                z = row['z']
                lay = row['lay']
                vertex_indices = row['vertex_indices']
                m = len(vertex_indices)

                # Create a space-separated string of vertex
                # indices
                vertex_indices_str = " ".join(
                    map(str, vertex_indices))

                f.write(f"{inode} {x} {y} {z} {lay} {m} "
                        f"{vertex_indices_str}\n")

        print(
            f"Successfully wrote grid data to {filename}")

    except IOError as e:
        print(f"Error writing to file: {e}")


if __name__ == "__main__":
    vertices_df, node_data_df = read_gsf_file(r"C:\workspace\codes\gridutil\tests\simple_tutorial_voronoi.gsf")

    write_gsf_file(vertices_df,node_data_df, filename="voronoi_grid.gsf")


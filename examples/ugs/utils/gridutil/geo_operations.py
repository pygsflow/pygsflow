import pandas as pd
import shapely.geometry as geom

def compute_buffer(stream_lines, buffer_distance):
    """Compute buffer around stream lines and merge into one polygon."""
    from shapely.ops import unary_union
    
    # Compute buffer around stream lines
    stream_buffers = []

    if isinstance(stream_lines, list):
        # If stream_lines is a list of LineStrings
        for line in stream_lines:
            buffer_polygon = line.buffer(buffer_distance)
            stream_buffers.append(buffer_polygon)
    else:
        # If stream_lines is a single LineString
        buffer_polygon = stream_lines.buffer(buffer_distance)
        stream_buffers.append(buffer_polygon)

    print(f"Created {len(stream_buffers)} buffer polygons around stream lines")
    
    # Merge all buffer polygons into one valid polygon
    if len(stream_buffers) > 1:
        merged_buffer = unary_union(stream_buffers)
        print(f"Merged {len(stream_buffers)} polygons into one with area: {merged_buffer.area:.2f}")
    else:
        merged_buffer = stream_buffers[0]
        print(f"Single buffer polygon with area: {merged_buffer.area:.2f}")
    
    return merged_buffer

def intersect_lines_polygons(polygons, lines):
    """Find which lines intersect with polygons using spatial index for speed."""
    from shapely.strtree import STRtree
    
    intersecting_lines = []
    intersecting_polys = []
    
    # Create a mapping from polygon objects to their indices
    polygon_to_index = {id(poly): i for i, poly in enumerate(polygons)}
    
    # Build spatial index for polygons
    polygon_tree = STRtree(polygons)
    
    for i, line in enumerate(lines):
        # Find potential intersecting polygons using spatial index
        potential_polygons = polygon_tree.query(line)
        
        # Check actual intersection with potential candidates
        for polygon in potential_polygons:
            if line.intersects(polygon):
                poly_index = polygon_to_index[id(polygon)]
                intersecting_lines.append(i)
                intersecting_polys.append(poly_index)
                break  # Found intersection, move to next line
    
    df = pd.DataFrame(columns=['line_index', 'poly_index'])
    df['line_index'] = intersecting_lines
    df['poly_index'] = intersecting_polys
    return df

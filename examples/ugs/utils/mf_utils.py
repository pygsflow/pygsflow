import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from shapely.geometry import LineString, Point, Polygon
from shapely.strtree import STRtree

from utils.gridutil import plot_grid

DEBUG = False


def get_sfr_lines(mf):
    """Extract SFR reach-to-reach line segments from a model.

    Connects the centroids of consecutive reaches within each
    stream segment to form LineString geometries.

    Args:
        mf: Flopy MODFLOW model object with an SFR package.

    Returns:
        DataFrame with columns: segment, ireach,
        next_ireach, lineGeometry.
    """
    reach_data = mf.sfr.reach_data
    df_reach_data = pd.DataFrame(reach_data)

    cell_size = mf.modelgrid.delr[0]
    line_info = []

    for segment in df_reach_data['iseg'].unique():
        df_segment = df_reach_data[
            df_reach_data['iseg'] == segment
        ].sort_values('ireach')

        # Connect consecutive reach centroids
        for i in range(len(df_segment) - 1):
            cur = df_segment.iloc[i]
            nxt = df_segment.iloc[i + 1]

            cur_center = (
                (cur['j'] - 1) * cell_size + cell_size / 2,
                (cur['i'] - 1) * cell_size + cell_size / 2
            )
            nxt_center = (
                (nxt['j'] - 1) * cell_size + cell_size / 2,
                (nxt['i'] - 1) * cell_size + cell_size / 2
            )

            line = LineString([cur_center, nxt_center])
            line_info.append([
                segment, cur['ireach'],
                nxt['ireach'], line
            ])

    lines_df = pd.DataFrame(
        line_info,
        columns=[
            'segment', 'ireach',
            'next_ireach', 'lineGeometry'
        ]
    )
    return lines_df


def get_unique_vertices(dis):
    """Find unique vertices in a structured MODFLOW grid.

    Iterates over all grid corners and collects the unique
    (x, y, z) coordinates from the top and bottom arrays.
    Vertices at interior corners shared by multiple cells
    are stored only once.

    Args:
        dis: Flopy discretization object with delr, delc,
            nrow, ncol, nlay, top, and botm attributes.

    Returns:
        numpy.ndarray of shape (n_vertices, 3) with x, y, z
        coordinates.
    """
    delr = dis.delr.array
    delc = dis.delc.array
    nrow = dis.nrow
    ncol = dis.ncol
    nlay = dis.nlay
    top = dis.top.array
    bottom = dis.botm.array

    # Edge coordinates along each axis
    x_coords = np.concatenate([[0], np.cumsum(delr)])
    y_coords = np.concatenate([[0], np.cumsum(delc)])

    # Pre-allocate for worst-case vertex count
    max_vertices = (nrow + 1) * (ncol + 1) * (nlay + 2)
    vertices_x = np.empty(max_vertices, dtype=np.float64)
    vertices_y = np.empty(max_vertices, dtype=np.float64)
    vertices_z = np.empty(max_vertices, dtype=np.float64)
    vertex_count = 0

    for i in range(nrow + 1):
        for j in range(ncol + 1):
            x = x_coords[j]
            y = y_coords[i]

            z_values = []
            if i < nrow and j < ncol:
                z_values.append(top[i, j])
                z_values.extend(bottom[:, i, j])

            # Keep only unique elevations at this (x, y)
            z_values = np.unique(z_values)

            n_z = len(z_values)
            vertices_x[vertex_count:vertex_count + n_z] = x
            vertices_y[vertex_count:vertex_count + n_z] = y
            vertices_z[vertex_count:vertex_count + n_z] = z_values
            vertex_count += n_z

    vertices = np.column_stack([
        vertices_x[:vertex_count],
        vertices_y[:vertex_count],
        vertices_z[:vertex_count]
    ])
    return vertices


def intersect_line_with_grid(mi):
    """Intersect SFR stream lines with the unstructured grid.

    Uses a spatial index (STRtree) to efficiently find which
    grid cells each stream segment passes through. Records
    the intersection length and distance from the segment's
    upstream end for reach ordering.

    Results are stored in mi.usg_sfr_df with columns:
    segment_id, node_id, distance_from_segment_upstream,
    intersection_length.

    Args:
        mi: Model_info object with sfr_df (containing
            seg_line column) and oct_grid2d (list of
            Shapely polygons).
    """
    lines = mi.sfr_df["seg_line"].values
    polygons = mi.oct_grid2d
    tree = STRtree(polygons)

    # Minimum intersection length to count as a valid reach
    min_length_threshold = 5

    usg_sfr = []
    sfr_polygons = []

    for i, line in enumerate(lines):
        candidates = tree.query(line, predicate='intersects')
        start_point = Point(line.xy[0][0], line.xy[1][0])

        for poly_idx in candidates:
            intersection = line.intersection(polygons[poly_idx])

            if intersection.length >= min_length_threshold:
                sfr_polygons.append(polygons[poly_idx])
                centroid = polygons[poly_idx].centroid
                dist = start_point.distance(centroid)
                usg_sfr.append([
                    i + 1, poly_idx,
                    dist, intersection.length
                ])

    mi.usg_sfr_df = pd.DataFrame(
        usg_sfr,
        columns=[
            'segment_id', 'node_id',
            'distance_from_segment_upstream',
            'intersection_length'
        ]
    )

    if DEBUG:
        fig, ax = plt.subplots()
        plot_grid(polygons, ax=ax, edge_color='green')
        plot_grid(sfr_polygons, ax=ax, edge_color='red')
        for line in lines:
            ax.plot(line.xy[0], line.xy[1], 'k-')
        plt.show()

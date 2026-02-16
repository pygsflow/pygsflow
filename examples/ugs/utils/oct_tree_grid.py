import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from shapely.geometry import Polygon, LineString
from shapely.ops import unary_union

import flopy
from flopy.utils import flopy_io
from flopy.utils.gridgen import Gridgen

from utils import gridutil, prms_utils


def assert_gridgen_executable():
    """Check that the gridgen executable is on PATH."""
    gridgen_exe = flopy.which("gridgen")
    if gridgen_exe is None:
        print(
            "Warning, gridgen is not in your path. "
            "You will need to provide a full path to "
            "the gridgen binary executable."
        )
    else:
        print(
            f"gridgen executable was found at: "
            f"{flopy_io.relpath_safe(gridgen_exe)}"
        )


def get_line_coords(sfr):
    """Extract stream segment LineStrings from SFR data.

    Builds a LineString for each segment by connecting reach
    cell centres in order. Downstream segments are connected
    by appending the first point of the outseg.

    Args:
        sfr: Flopy SFR package object.

    Returns:
        DataFrame with columns: iseg, seg_line (LineString),
        outseg, ij_coords.
    """
    nrow = sfr.parent.nrow
    ncol = sfr.parent.ncol
    cell_size = sfr.parent.modelgrid.delr[0]
    delr = sfr.parent.modelgrid.delr
    delc = sfr.parent.modelgrid.delc
    segment_data = pd.DataFrame(sfr.segment_data[0])
    reach_data = pd.DataFrame(sfr.reach_data)
    unique_iseg = segment_data["nseg"].unique()

    reach_info = []
    for iseg in unique_iseg:
        if iseg <= 0:
            continue
        cur = reach_data[reach_data["iseg"] == iseg]
        reach_ids = np.sort(cur["ireach"].values)

        seg_line = []
        ij_coords = []
        for ireach in reach_ids:
            row = cur[cur["ireach"] == ireach]
            i = row["i"].values
            j = row["j"].values
            # Cell centre with small offset to avoid edge issues
            x1 = delr[j] * j - delr[j] / 2
            y1 = delc.sum() - (delc[i] * i - delc[i] / 2)
            eps = cell_size / 10000
            seg_line.append((x1 + eps, y1 + eps))
            ij_coords.append((i[0], j[0]))

        outseg = segment_data[
            segment_data["nseg"] == iseg
        ]["outseg"].values[0]

        reach_info.append(
            [iseg, seg_line, outseg, ij_coords]
        )

    df = pd.DataFrame(
        reach_info,
        columns=["iseg", "seg_line", "outseg", "ij_coords"]
    )

    # Connect segment endpoints to downstream segment starts
    for iseg in unique_iseg:
        df_seg = df[df["iseg"] == iseg]
        seg_line = df_seg["seg_line"].values[0]
        outseg = df_seg["outseg"].values[0]

        if outseg > 0:
            df_outseg = df[df["iseg"] == outseg]
            first_point = df_outseg["seg_line"].values[0][0]
            seg_line.append(first_point)

        df.loc[
            df["iseg"] == iseg, "seg_line"
        ] = LineString(seg_line)

    return df


def mask_to_polygon(mask, modelgrid):
    """Convert a boolean watershed mask to a Gridgen polygon.

    Uses matplotlib contour tracing to extract the boundary
    of the mask, scaled to model coordinates.

    Args:
        mask: 2D boolean array (True = watershed).
        modelgrid: Flopy StructuredGrid with delr/delc.

    Returns:
        Polygon in Gridgen format: [[[vertices]]], or []
        if no contour is found.
    """
    fig, ax = plt.subplots()

    nrow, ncol = mask.shape
    row_numbers = np.arange(nrow)
    col_numbers = np.arange(ncol)

    row_grid, col_grid = np.meshgrid(
        row_numbers, col_numbers, indexing='ij'
    )
    row_grid = row_grid * modelgrid.delr[0]
    col_grid = col_grid * modelgrid.delc[0]

    contours = ax.contour(
        col_grid, np.flipud(row_grid),
        mask.astype(float), levels=[0.5]
    )

    polygon_coords = []
    if hasattr(contours, 'allsegs') and contours.allsegs:
        all_paths = contours.allsegs[0]
        if all_paths:
            largest_path = max(
                all_paths, key=lambda p: len(p)
            )
            polygon_coords = [
                tuple(v) for v in largest_path.tolist()
            ]

    plt.close(fig)

    if polygon_coords:
        return [[polygon_coords]]

    print("Warning: No contour found for watershed mask")
    return []


def create_buffer_polygons_around_streams(
    stream_lines, buffer_distance, cell_size,
    active_domain_polygon=None, inset_distance=10.0
):
    """Create dissolved buffer polygons around stream lines.

    Buffers each stream line and merges them into a single
    polygon. If an active domain is provided, the buffer is
    clipped to stay within an inset boundary.

    Args:
        stream_lines: List of Shapely LineString objects.
        buffer_distance: Buffer width in model units.
        cell_size: Grid cell size (for fallback method).
        active_domain_polygon: Gridgen-format polygon to
            clip against (optional).
        inset_distance: Distance to shrink the active domain
            boundary before clipping (default 10.0).

    Returns:
        Buffer polygon(s) in Gridgen format: [[[vertices]]].
    """
    individual_buffers = []
    for line in stream_lines:
        individual_buffers.append(line.buffer(buffer_distance))

    if not individual_buffers:
        return []

    dissolved = unary_union(individual_buffers)

    # Clip to active domain if provided
    if active_domain_polygon and len(active_domain_polygon) > 0:
        domain_coords = active_domain_polygon[0][0]
        domain_poly = Polygon(domain_coords)
        inset_poly = domain_poly.buffer(-inset_distance)

        if not inset_poly.is_empty:
            clipped = dissolved.intersection(inset_poly)
            if not clipped.is_empty:
                dissolved = clipped

    # Convert to Gridgen format
    if hasattr(dissolved, 'exterior'):
        return [[list(dissolved.exterior.coords)]]
    else:
        result = []
        for geom in dissolved.geoms:
            result.append([list(geom.exterior.coords)])
        return result


def grid_2d_from_gridgen(mi):
    """Extract 2D cell polygons from the Gridgen output.

    Reads the unstructured grid vertices and cell
    connectivity, then constructs a Shapely Polygon for
    each cell.

    Args:
        mi: Model_info object with gridgen_obj attribute.

    Returns:
        List of Shapely Polygon objects (one per cell).
    """
    g = mi.gridgen_obj
    mi.gridprops = g.get_gridprops_disu5()
    mi.unstructured_grid = g.get_gridprops_unstructuredgrid()

    vertices = np.array(mi.unstructured_grid['vertices'])
    grid2d = []
    for cell in mi.unstructured_grid['iverts']:
        ver_idx = cell[:-1]
        poly = Polygon(
            [(x, y) for x, y in vertices[ver_idx, 1:]]
        )
        grid2d.append(poly)
    return grid2d


def grid_from_structured_dis(dis):
    """Build a cell DataFrame from a structured DIS package.

    Creates a row per cell with node_id, Shapely polygon,
    top/bottom elevations, and layer index.

    Args:
        dis: Flopy discretization object.

    Returns:
        DataFrame with columns: node_id, geometry, top,
        botm, layer.
    """
    cols = ["node_id", "geometry", "top", "botm", "layer"]
    rows = []
    modelgrid = dis.parent.modelgrid

    for k in range(dis.nlay):
        for i in range(dis.nrow):
            for j in range(dis.ncol):
                node_id = modelgrid.get_node((k, i, j))
                verts = modelgrid.get_cell_vertices(node_id)
                rows.append({
                    "node_id": node_id[0],
                    "geometry": Polygon(verts),
                    "top": modelgrid.top[i, j],
                    "botm": modelgrid.botm[k, i, j],
                    "layer": k,
                })

    return pd.DataFrame(rows, columns=cols)


def _build_gsf_data(mi, vertices):
    """Build GSF vertex and node DataFrames from Gridgen output.

    Creates 3D vertices by duplicating 2D vertices at the top
    and bottom elevations, and builds node data with centroid
    coordinates and vertex index lists.

    Args:
        mi: Model_info object with unstructured_grid and
            gridprops attributes.
        vertices: DataFrame of 2D vertices (id, x, y).

    Returns:
        Tuple of (vertices_3d DataFrame, node_data DataFrame).
    """
    n_verts_top = vertices.shape[0]

    # Create top and bottom vertex copies
    verts_top = pd.DataFrame(
        np.array(mi.unstructured_grid['vertices']),
        columns=['id', 'x', 'y']
    )
    verts_top['z'] = None
    verts_top['id'] = verts_top['id'].astype(int)
    verts_top.sort_values(by='id', inplace=True)

    verts_bot = verts_top.copy()

    # Assign z-values from cell top/bottom elevations
    iverts2 = []
    for jj, iv in enumerate(mi.unstructured_grid['iverts']):
        cell_verts = iv[:-1]
        z_top = mi.unstructured_grid['top'][jj]
        z_bot = mi.unstructured_grid['botm'][jj]

        verts_top.loc[cell_verts, 'z'] = z_top
        verts_bot.loc[cell_verts, 'z'] = z_bot

        # Bottom vertex IDs are offset by n_verts_top
        ids_bot = [v + n_verts_top for v in cell_verts]
        # Convert to 1-based indices
        new_iverts = [v + 1 for v in cell_verts + ids_bot]
        iverts2.append(new_iverts)

    vertices_3d = pd.concat([verts_top, verts_bot])
    vertices_3d = pd.DataFrame(
        vertices_3d, columns=['x', 'y', 'z']
    )

    # Build node data
    node_data = pd.DataFrame(
        columns=[
            'inode', 'x', 'y', 'z', 'lay',
            'vertex_indices'
        ]
    )
    ug = mi.unstructured_grid
    node_data['inode'] = 1 + np.arange(mi.gridprops['nodes'])
    node_data['x'] = ug['xcenters']
    node_data['y'] = ug['ycenters']
    node_data['z'] = (ug['top'] + ug['botm']) / 2.0
    node_data['lay'] = 1
    node_data['vertex_indices'] = iverts2

    return vertices_3d, node_data


def create_oct_tree_grid(mi, stream_buffer_distance=100.0):
    """Generate an octree-refined unstructured grid.

    Full workflow:
    1. Verify gridgen executable
    2. Extract watershed boundary and stream lines
    3. Set active domain and add stream-buffer refinement
    4. Run Gridgen to build the octree grid
    5. Compute structured-to-unstructured cell mapping
    6. Write the grid specification (GSF) file

    Args:
        mi: Model_info object with fine_gsf, gridgen_exe,
            and path attributes.
        stream_buffer_distance: Buffer width around streams
            for grid refinement (model units, default 100).
    """
    assert_gridgen_executable()

    gridgen_ws = os.path.join(mi.ws, "gridgen")
    os.makedirs(gridgen_ws, exist_ok=True)
    print(
        f"Model workspace is : "
        f"{flopy_io.scrub_login(mi.ws)}"
    )
    print(
        f"Gridgen workspace is : "
        f"{flopy_io.scrub_login(gridgen_ws)}"
    )

    ms = mi.fine_gsf.mf
    dis = ms.dis
    source_grid = grid_from_structured_dis(dis)

    g = Gridgen(
        ms.modelgrid,
        model_ws=gridgen_ws,
        exe_name=mi.gridgen_exe
    )

    # Define active domain from watershed mask
    mi.watershed = prms_utils.get_prms_watershed(mi.fine_gsf)
    polygon = mask_to_polygon(
        mi.watershed, ms.modelgrid
    )
    g.add_active_domain(polygon, range(mi.nlay))

    # Extract stream line geometries
    sfr = ms.sfr
    sfr_df = get_line_coords(sfr)
    mi.sfr_df = sfr_df
    stream_lines = sfr_df["seg_line"].values.tolist()

    # Add stream-buffer refinement zone
    cell_size = ms.modelgrid.delr[0]
    inset_distance = 10.0
    buffer_polygons = create_buffer_polygons_around_streams(
        stream_lines,
        stream_buffer_distance,
        cell_size,
        polygon,
        inset_distance
    )

    if buffer_polygons:
        print(
            f"Created {len(buffer_polygons)} buffer "
            f"polygons around streams"
        )
        g.add_refinement_features(
            buffer_polygons, "polygon", 2, range(mi.nlay)
        )

    # Build the octree grid
    g.build(verbose=True)

    mi.gridgen_obj = g
    mi.oct_grid2d = grid_2d_from_gridgen(mi)

    # Compute structured-to-unstructured area mapping
    struct_geom = source_grid[
        source_grid['layer'] == 0
    ]['geometry'].values.tolist()
    mi.mapping_df = gridutil.compute_grid_exact_intersections(
        source_grid=struct_geom,
        target_grid=mi.oct_grid2d
    )

    # Build and write GSF file
    vertices = pd.DataFrame(
        np.array(mi.unstructured_grid['vertices']),
        columns=['id', 'x', 'y']
    )
    vertices_3d, node_data = _build_gsf_data(mi, vertices)

    gsf_path = os.path.join(
        mi.usg_model_ws, "grid_spec.gsf"
    )
    gridutil.write_gsf_file(
        vertices_3d, node_data, filename=gsf_path
    )

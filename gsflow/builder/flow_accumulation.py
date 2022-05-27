import numpy as np
from collections import OrderedDict
import random
from ..utils.sfr_renumber import Topology
from ..utils import gsflow_io
from . import Defaults
from flopy.modflow import ModflowSfr2
import shapefile
import inspect


class FlowAccumulation(object):
    """
    Class to perform flow accumulation on a DEM (raster) or resampled to
    a model grid

    Parameters
    ----------
    data : np.ndarray
        two dimensional numpy array of dem data
    xcenters : np.ndarray
        a two dimensional array of x coordinate cell centers
    ycenters : np.ndarray
        two dimensional array of y coordinate cell centers
    acc_type : str
        flow accumulation type, currently supported options are "d8"
    hru_type : np.array
        optional numpy array of hru type numbers that can be used as
        a masking array to exclude lakes, swales, etc...
        0 == inactive
        1 == land (included in flow accumulation)
        2 == lake (excluded from flow accumulation)
        3 == swale (excluded from flow accumulation)
    closed_basin : bool
        method to indicate that basin is closed without considering lake
        flow. If true hru_type 2 is used in the flow direction calculations.
        False ignores hru_type 2.
    flow_dir_array : np.ndarray
        previously calculated flow direction array of dimension nrow, ncol
        that can be supplied to FlowAccumulation
    verbose : bool
        flag to print verbose output

    """

    def __init__(
        self,
        data,
        xcenters,
        ycenters,
        acc_type="d8",
        hru_type=None,
        closed_basin=False,
        flow_dir_array=None,
        verbose=False,
    ):

        self._defaults = Defaults().to_dict()

        # set the fa type
        self._acc_type = acc_type
        # flow directions vectors set up as top to bot, left to right
        self._d8_vectors = np.array([32, 64, 128, 16, 1, 8, 4, 2], dtype=int)
        self._quiver_u = {
            32: -1,
            64: 0,
            128: 1,
            16: -1,
            1: 1,
            8: -1,
            4: 0,
            2: 1,
            -1: np.nan,
            -2: np.nan,
        }
        self._quiver_v = {
            32: 1,
            64: 1,
            128: 1,
            16: 0,
            1: 0,
            8: -1,
            4: -1,
            2: -1,
            -1: np.nan,
            -2: np.nan,
        }
        self._inner_idx = []
        self._graph = OrderedDict()
        self._visited = []
        self._stack = []
        self._dest = None
        self._solved_flats = {}

        # adds an extra row and col of data
        data, xcenters, ycenters, hru_type = self._buffer_data(
            data, xcenters, ycenters, hru_type
        )

        self._offset = data.shape[1]
        self._offsets = np.array(
            [
                (-1 * self._offset) - 1,
                (-1 * self._offset),
                (-1 * self._offset) + 1,
                -1,
                1,
                self._offset - 1,
                self._offset,
                self._offset + 1,
            ]
        )
        self._shape = data.shape
        self._data = data.ravel()
        self._xcenters = xcenters.ravel()
        self._ycenters = ycenters.ravel()
        self._hru_type = hru_type.ravel()
        self._flow_directions = np.ones((self._data.size,)) * -1
        self._closed_basin = closed_basin

        if flow_dir_array is None:

            self._flow_directions = np.ones((self._data.size,)) * -1
            self._closed_basin = closed_basin
            if self._closed_basin:
                self._cb_inner_idx = []
                self._cb_flow_directions = np.ones((self._data.size,)) * -1

        else:
            if self._closed_basin:
                err = (
                    "Supplying flow direction array not supported "
                    "for closed basin."
                )
                raise Exception(err)

            self._flow_directions = np.pad(
                flow_dir_array, 1, "constant", constant_values=-1
            ).ravel()
            self._set_inners()

        # for flow accumulation
        self._dir_coords = {
            32: (-1, -1),
            64: (-1, 0),
            128: (-1, 1),
            16: (0, -1),
            1: (0, 1),
            8: (1, -1),
            4: (1, 0),
            2: (1, 1),
        }

        # for watersheds
        self._d8_vectors_r = np.array(
            list(reversed(self._d8_vectors)), dtype=int
        )

        self._offset_dict = {
            32: self._offsets[0],
            64: self._offsets[1],
            128: self._offsets[2],
            16: self._offsets[3],
            1: self._offsets[4],
            8: self._offsets[5],
            4: self._offsets[6],
            2: self._offsets[7],
        }

        self._size = self._data.size
        self._dijkstra = False
        self.verbose = verbose
        self._wpp = None

    def get_dem_data(self):
        """
        Method to get processed DEM data after flow accumulation

        Returns
        -------
            np.ndarray
        """
        return np.reshape(self._data, self._shape)[1:-1, 1:-1]

    def get_hru_type(self):
        return np.reshape(self._hru_type, self._shape)[1:-1, 1:-1]

    def flow_directions(self, dijkstra=False, breach=0.0):
        """
        Method to get flow directions array

        Parameters
        ----------
        dijkstra : bool
            method to use a modified dijkstra algorithmn to solve
            tricky flat areas. Default is False and uses a distance based
            topology method to solve flat areas.
        breach : float
            absolute value of breaching tolerance for digital dams. Use
            caution while applying breaching values. These should be small
            numbers.

        Returns
        -------
            np.ndarray of flow directions
        """
        self._dijkstra = dijkstra
        if np.all(self._flow_directions == -1):
            # builds inner indices
            if self._closed_basin:
                self._set_inners(cb=True)
            else:
                self._set_inners()

            # fills pits in DEM
            self._fill_pits()

            self._d8_flow_directions(breach)

            # reset the inners after performing the closed basin calculation
            if self._closed_basin:
                self._cb_inner_idx = list(np.copy(self._inner_idx))
                self._inner_idx = []
                self._set_inners()

        if self._closed_basin:
            self._cb_flow_directions = np.copy(self._flow_directions)
            self._flow_directions[self._hru_type == 2] = -1

        flow_directions = np.copy(self._flow_directions)
        flow_directions.shape = self._shape

        return flow_directions[1:-1, 1:-1]

    @property
    def get_vectors(self):
        """
        Method to get flow vectors array

        Returns
        -------
             (u, v): tuple (np.ndarray of flow vectors)
        """
        flow_direction = self.flow_directions()
        u = np.zeros(flow_direction.shape)
        v = np.zeros(flow_direction.shape)
        for i, row in enumerate(flow_direction):
            for j, val in enumerate(row):
                u[i, j] = self._quiver_u[val]
                v[i, j] = self._quiver_v[val]
        return u, v

    def _buffer_data(self, data, xcenters, ycenters, hru_type):
        """
        Method to prep data for flow direction calculations

        Parameters
        ----------
        data : np.ndarray
        xcenters : np.ndarray
        ycenters : np.ndarray
        hru_type : np.ndarray

        Returns
        -------
        buffered arrays (data, xcenters, ycenters, hru_type)

        """
        shape = data.shape
        xsize = shape[0] + 2
        ysize = shape[1] + 2

        new_xcenters = np.zeros((xsize, ysize))
        new_ycenters = np.zeros((xsize, ysize))
        new_data = np.zeros((xsize, ysize))
        new_hru_type = np.ones((xsize, ysize), dtype=int)

        xdist = np.abs(xcenters[0, 0] - xcenters[0, 1])
        ydist = np.abs(ycenters[0, 0] - ycenters[1, 0])
        fill_val = np.max(data) + 1e05

        # buffer the x-centers array
        new_xcenters[1:-1, 1:-1] = xcenters
        new_xcenters[0, :] = new_xcenters[1, :]
        new_xcenters[:, 0] = new_xcenters[:, 1] - xdist
        new_xcenters[-1, :] = new_xcenters[-2, :]
        new_xcenters[:, -1] = new_xcenters[:, -2] + xdist

        # buffer the y-centers array
        new_ycenters[1:-1, 1:-1] = ycenters
        new_ycenters[0, :] = new_ycenters[1, :] + ydist
        new_ycenters[:, 0] = new_ycenters[:, 1]
        new_ycenters[-1, :] = new_ycenters[-2, :] - ydist
        new_ycenters[:, -1] = new_ycenters[:, -2]

        new_data[:, :] = fill_val
        new_data[1:-1, 1:-1] = data

        if hru_type is not None:
            new_hru_type[1:-1, 1:-1] = hru_type

        return new_data, new_xcenters, new_ycenters, new_hru_type

    def _set_inners(self, cb=False):
        """
        Method to set inner indicies for flow direction
        and accumulation calculations.

        """
        calc_hru = (1,)
        if cb:
            calc_hru = (1, 2)

        nnodes = self._shape[0] * self._shape[1]
        for ix in range(nnodes):
            if ix < self._offset:
                continue
            elif ix % self._offset == 0:
                continue
            elif (ix + 1) % self._offset == 0:
                continue
            elif ix > (nnodes - self._offset):
                continue
            elif self._hru_type[ix] not in calc_hru:
                continue
            else:
                self._inner_idx.append(ix)

    def _get_fa_stack_node(self):
        """
        Method to pop a node of the stack

        Returns
        -------
            int : node number
        """
        if len(self._stack) > 0:
            return self._stack.pop(0)
        else:
            return None

    def _check_if_solved(self, dir_idx, node_numbers):
        """
        Method to check if this flat area has a solution

        Parameters
        ----------
        dir_idx : list
            flow direction index
        node_numbers : list
            list of node numbers

        Returns
        -------
            dir_idx, node_numbers
        """
        dist = 999
        solved = None
        s_ix = None
        for ix, node in enumerate(node_numbers):
            if node in self._solved_flats:
                tdist = self._solved_flats[node]
                if tdist < dist:
                    solved = node
                    s_ix = ix
        if solved is not None:
            node_numbers = np.array([solved], dtype=int)
            dir_idx = np.array([dir_idx[s_ix]], dtype=int)

        return dir_idx, node_numbers

    def _resolve_flats(self, ix0, dem, breach=0.0):
        """
        Non recursive method that builds a topology and solves
        based on imporved distance of connectivity

        Parameters
        ----------
        ix0 : int
            initial index/node number
        dem : np.ndarray
            digital elevation map to calculate drops

        Returns
        -------
            int : flow direction vector

        """
        calc_hru = (1,)
        if self._closed_basin:
            calc_hru = (1, 2)

        while self._stack:
            ix = self._get_fa_stack_node()
            if self._hru_type[ix] not in calc_hru:
                continue

            if ix in self._graph:
                # we have visited this node
                continue

            dir_idx, node_numbers = self._calculate_drop(ix, dem, breach)
            if len(dir_idx) == 1:
                tdem = dem[node_numbers[0]]
                if self._dest is None:
                    self._dest = [node_numbers[0], tdem, 1]
                else:
                    if self._dest[-1] == 1:
                        if tdem < self._dest[1]:
                            self._dest = [node_numbers[0], tdem, 1]
                    else:
                        self._dest = [node_numbers[0], tdem, 1]
            else:
                diff = dir_idx[1:] - dir_idx[:-1]
                diff = diff == 1
                consecutive = False
                if np.alltrue(diff):
                    consecutive = True

                if not self._dijkstra:
                    if len(dir_idx) == 2 and consecutive:
                        if self._dest is None:
                            ii = random.randint(0, 1)
                            tdem = dem[node_numbers[ii]]
                            self._dest = [node_numbers[ii], tdem, 2]
                        elif self._dest[-1] > 2:
                            ii = random.randint(0, 1)
                            tdem = dem[node_numbers[ii]]
                            self._dest = [node_numbers[ii], tdem, 2]

                    elif len(dir_idx) == 3 and consecutive:
                        if self._dest is None:
                            ii = random.randint(0, 2)
                            tdem = dem[node_numbers[ii]]
                            self._dest = [node_numbers[1], tdem, 3]
                        elif self._dest[-1] > 3:
                            ii = random.randint(0, 2)
                            tdem = dem[node_numbers[ii]]
                            self._dest = [node_numbers[1], tdem, 3]

            self._graph[ix] = node_numbers
            self._stack += list(node_numbers)

        if self._dest is None:
            # we should then set all of these nodes to -2
            for node in self._graph.keys():
                self._flow_directions[node] = -2
            return -2

        if self._dijkstra:
            print("Applying Dijkstra solution to resolve flat area:")
            stack = [self._dest[0]]
            destination = self._dest[0]
            trace_outlet = {}
            n = 0
            dn = 0
            graph_nodes = list(self._graph.keys())
            while len(trace_outlet) < len(graph_nodes):
                if self.verbose:
                    if n % 20 == 0:
                        print(f"pass # {n}")
                if n > len(graph_nodes):
                    # this seems to work even when n > len(graph_nodes)
                    # consider removing unless problems arise
                    pass

                tdest = []
                for ix, node_numbers in self._graph.items():
                    for dest in stack:
                        if dest in node_numbers:
                            if ix not in trace_outlet:
                                if dest == destination:
                                    dist = 0
                                else:
                                    if dest not in trace_outlet:
                                        dist = 999
                                    else:
                                        dist = trace_outlet[dest][1] + 1
                                trace_outlet[ix] = [dest, dist]
                                tdest.append(ix)
                            else:
                                dist = trace_outlet[dest][1] + 1
                                prev_dist = trace_outlet[ix][1]
                                if dist < prev_dist:
                                    trace_outlet[ix] = [dest, dist]

                for dest in tdest:
                    if dest in self._graph:
                        self._graph.pop(dest)

                if not tdest:
                    # goto the next node number
                    dn += 1
                    try:
                        node = graph_nodes[dn]
                    except IndexError:
                        print("break")
                    stack = [node]
                else:
                    stack = tdest
                n += 1

            rtn_val = -2
            for ix, node_number in trace_outlet.items():
                node = node_number[0]
                idxs = self._offsets + ix
                dir_idx = np.where(idxs == node)[0]
                vector = self._d8_vectors[dir_idx[0]]
                if ix == ix0:
                    rtn_val = vector
                self._flow_directions[ix] = vector
                self._solved_flats[ix] = node_number[1]

        else:
            rtn_val = -2
            for ix, node_numbers in sorted(self._graph.items()):
                # use distance equation to determine the proper direction
                connection = self._minimize_distance(node_numbers)
                idxs = self._offsets + ix
                dir_idx = np.where(idxs == connection)[0]
                vector = self._d8_vectors[dir_idx[0]]
                if ix == ix0:
                    rtn_val = vector
                self._flow_directions[ix] = vector

        return rtn_val

    def _minimize_distance(self, node_numbers):
        """
        Method to implement distance oriented graph
        to find most direct route from outlet to inlet
        Inspired by Dijkstra's algorithm

        Parameters
        ----------
        node_numbers : list or np.ndarray
            list of node numbers connected to a cell

        Returns
        -------
            int : connection, node number provides the most reduced distance to
               a destination cell

        """
        dxdy = (self._xcenters[self._dest[0]], self._ycenters[self._dest[0]])
        nx = np.array([self._xcenters[i] for i in node_numbers])
        ny = np.array([self._ycenters[i] for i in node_numbers])

        a2 = (nx - dxdy[0]) ** 2
        b2 = (ny - dxdy[1]) ** 2
        distance = np.sqrt(a2 + b2)

        ix = np.where(distance == np.min(distance))[0]
        connection = node_numbers[ix[0]]
        return connection

    def _d8_flow_directions(self, breach=0.0):
        """
        Method to create D8 flow directions

        Parameters
        ----------
        breach : float
            absolute threshold to ignore/breach digital dams
        """
        dem = self._data
        for ix in self._inner_idx:
            if self._flow_directions[ix] != -1:
                continue
            self._graph = OrderedDict()
            self._stack = []
            self._dest = None
            vector = self.__d8_flow_directions(ix, dem, breach)
            self._flow_directions[ix] = vector

    def __d8_flow_directions(self, ix, dem, breach=0.0):
        """
        Recursive calculation method that allows for "conflict
        resolution", ex. flat areas in the basin or equal elevation
        differences.

        Parameters:
        ----------
        ix : int
            index of cell to set flow direction
        offsets : np.array
            set of array offsets to find nearest neighbors
        breach : float
            absolute threshold to ignore/breah digital dams

        Returns
        -------
            int : flow direction
        """
        dir_idx, node_numbers = self._calculate_drop(ix, dem, breach)

        if dir_idx.size == 1:
            ii = dir_idx[0]
            return self._d8_vectors[ii]
        else:
            dir_idx = np.sort(dir_idx)
            diff = dir_idx[1:] - dir_idx[:-1]
            diff = diff == 1
            consecutive = False
            if np.alltrue(diff):
                consecutive = True

            if not self._dijkstra:
                if dir_idx.size == 2 and consecutive:
                    ii = random.randint(0, 1)
                    return self._d8_vectors[dir_idx[ii]]

                else:
                    if dir_idx.size == 3:
                        if consecutive:
                            return self._d8_vectors[dir_idx[1]]

            else:
                pass

            dir_idx, node_numbers = self._check_if_solved(
                dir_idx, node_numbers
            )

            if len(dir_idx) == 1:
                ii = dir_idx[0]
                return self._d8_vectors[ii]

            self._graph[ix] = node_numbers
            self._stack += list(node_numbers)
            vector = self._resolve_flats(ix, dem, breach)
            return vector

    def _calculate_drop(self, ix, dem, breach=0.0):
        """
        Method to calculate the drop and return the indicies
        and node_numbers of the greatest drop

        ix : int
            index number
        dem : np.ndarray
            array of elevations
        breach : float
            breaching value to ignore digital dams

        Returns
        -------
            dir_idx, node_numbers
        """
        idxs = self._offsets + ix

        cell_elevation = dem[ix]
        neighbor_elevation = dem[idxs]
        # prevent algortithm from considering inactive cells
        inactive = np.where(self._hru_type[idxs] == 0)[0]
        if len(inactive > 0):
            neighbor_elevation[inactive] = 1e20

        xcell = self._xcenters[ix]
        ycell = self._ycenters[ix]
        xneighbor = self._xcenters[idxs]
        yneighbor = self._ycenters[idxs]

        asq = (xneighbor - xcell) ** 2
        bsq = (yneighbor - ycell) ** 2
        dist = np.sqrt(asq + bsq)

        drop = neighbor_elevation - cell_elevation

        drop = np.where(np.abs(drop) <= breach, 0, drop)

        drop /= dist

        dir_idx = np.where(drop == np.min(drop))[0]
        node_numbers = idxs[dir_idx]

        return dir_idx, node_numbers

    def _fill_pits(self):
        """
        Method to fill pits in a DEM

        """
        dem = self._data
        offsets = self._offsets
        for ix in self._inner_idx:
            idxs = offsets + ix
            cell_elevation = dem[ix]
            neighbor_elvevation = dem[idxs]
            mask = neighbor_elvevation > cell_elevation
            if np.alltrue(mask):
                nval = np.min(neighbor_elvevation) + 1e-06
                dem[ix] = nval

        self._data = dem

    def flow_accumulation(self):
        """
        NIDP method to calculate the flow accumulation array using

        Returns
        -------
            np.ndarray of flow directions

        """

        # build NIDP array
        nidp_array = self._build_nidp()

        # set up flow directions
        flow_directions = self._flow_directions

        # flow accumulation array
        flow_accumulation = np.ones(nidp_array.shape)

        for ix in self._inner_idx:
            hru = self._hru_type[ix]
            if hru == 0:
                continue
            nidp_val = nidp_array[ix]
            if nidp_val != 0:
                continue
            n = ix
            naccu = 0
            while self._next_cell(n):
                flow_accumulation[n] += naccu
                naccu = flow_accumulation[n]
                if nidp_array[n] >= 2:
                    nidp_array[n] -= 1
                    break
                # get flow direction
                ndir = flow_directions[n]
                n = n + self._offset_dict[ndir]

        flow_accumulation.shape = self._shape
        return flow_accumulation[1:-1, 1:-1]

    def get_cascades(
        self, streams, pour_point=None, modelgrid=None, fmt="rowcol"
    ):
        """
        Method to calculate cascade parameters for PRMS

        Parameters
        ----------
        streams : _StreamsObj
            Stream information object returned from the
            FlowAccumulation.make_streams method
        pour_point : list, str
            optional, but required if watershed deliniation has not been run
            prior to running get_cascades. pour point can be supplied as either
            a shapefile, a list with an [(x, y)] tuple of coordinates,
            or as the model zero based [(row, column)] location.
        modelgrid : flopy.discretization.StructuredGrid

        fmt : str
            format of the pour point information. Acceptable types include
            "xy", "rowcol", and "shp"

        Returns
        -------
            _Cascades object
        """
        if pour_point is not None and self._wpp is None:
            if modelgrid is None:
                raise Exception(
                    "modelgrid must be supplied if pour_point is specified"
                )
            msg = (
                "running define watershed, watershed must be stored to "
                "run get_cascades"
            )

            gsflow_io._warning(
                msg, inspect.getframeinfo(inspect.currentframe())
            )
            self.define_watershed(pour_point, modelgrid, fmt=fmt)

        elif self._wpp is None:
            raise AssertionError("Watershed deliniation must be run first")

        hru_up_id, hru_down_id, hru_pct_up = self._build_nidp(cascades=True)

        # add the watershed outlet hru to cascades_ids
        hru_up_id.append(self._wpp)
        hru_down_id.append(self._wpp)
        hru_pct_up.append(1.0)

        # correct the hru_up_id and hru_down_id indicies to grid indicies
        hru_up_id = np.array(hru_up_id)
        hru_down_id = np.array(hru_down_id)

        hru_up_id = (
            (hru_up_id - self._shape[1])
            - (np.floor((hru_up_id - self._shape[1]) / self._shape[1]) * 2)
            - 1
        )
        hru_down_id = (
            (hru_down_id - self._shape[1])
            - (np.floor((hru_down_id - self._shape[1]) / self._shape[1]) * 2)
            - 1
        )
        hru_up_id = hru_up_id.astype(int)
        hru_down_id = hru_down_id.astype(int)
        hru_pct_up = np.array(hru_pct_up)

        iseg = streams.iseg.copy().ravel()
        strm_ix = []
        hru_strmseg_down_id = []
        for ix, hru_id in enumerate(hru_down_id):
            if iseg[hru_id] == 0:
                hru_strmseg_down_id.append(0)
            else:
                hru_strmseg_down_id.append(iseg[hru_id])
                strm_ix.append(ix)

        hru_strmseg_down_id = np.array(hru_strmseg_down_id)
        hru_up_id += 1
        hru_down_id += 1

        return _Cascades(
            hru_up_id, hru_down_id, hru_pct_up, hru_strmseg_down_id
        )

    def _next_cell(self, ix):
        """
        Method to check if code should advance to next cell in stack

        Parameters
        ----------
        ix : int

        Returns
        -------
            bool
        """
        nnodes = self._shape[0] * self._shape[1]
        if ix < self._offset:
            return False
        elif ix % self._offset == 0:
            return False
        elif (ix + 1) % self._offset == 0:
            return False
        elif ix > (nnodes - self._offset):
            return False
        elif self._hru_type[ix] != 1:
            return False
        elif self._flow_directions[ix] == -2:
            msg = (
                "sink present in data (fdir = -2), consider filling sinks "
                "in DEM and re-running flow direction and flow accumulation"
            )
            gsflow_io._warning(
                msg, inspect.getframeinfo(inspect.currentframe())
            )
            return False
        else:
            return True

    def _build_nidp(self, cascades=False):
        """
        Build the NIDP (Number of Incoming Dainage Paths)

        Parameters
        ----------
        cascades : bool
            if False this method calculates the NIDP, if True it calculates
            the hru_up_id and hru_down_id.

        Returns
        -------
            np.ndarray of NIDP or when cascades=True tuple of
            (hru_up_id, hru_down_id, hru_pct_up) arrays

        """
        nidp_array = np.zeros(self._flow_directions.shape)
        nidp_dir = [2, 4, 8, 1, 16, 128, 64, 32]
        # loop through inners
        hru_up_id = []
        hru_down_id = []
        hru_pct_up = []
        for ix in self._inner_idx:
            if self._hru_type[ix] == 0:
                continue
            # get neighbors
            if self._flow_directions[ix] == -1:
                raise Exception("Flow Direction not calculated")
            idxs = self._offsets + ix
            idxs_dir = self._flow_directions[idxs]

            if cascades:
                idxs_ix = np.where(nidp_dir == idxs_dir)[0]
                if len(idxs_ix) == 0:
                    continue
                else:
                    hru_up = list(idxs[idxs_ix])
                    # fix basin outlet issue for user supplied hru_type
                    if self._wpp in hru_up:
                        popix = hru_up.index(self._wpp)
                        hru_up.pop(popix)
                    hru_down = [ix] * len(hru_up)
                    hru_pct = [1.0 / len(hru_up)] * len(hru_up)

                hru_up_id += hru_up
                hru_down_id += hru_down
                hru_pct_up += hru_pct

            else:
                nidp = sum(
                    [
                        1
                        for nidpdir, idxdir in zip(nidp_dir, idxs_dir)
                        if nidpdir == idxdir
                    ]
                )
                nidp_array[ix] = nidp

        if cascades:
            return hru_up_id, hru_down_id, hru_pct_up
        else:
            return nidp_array

    def define_watershed(self, pour_point, modelgrid, fmt="rowcol"):
        """
        Method to perform watershed delineation on the flow direction array

        Parameters
        ----------
        pour_point : str, tuple
            pour point can be supplied as either a shapefile, a list with an
            [(x, y)] tuple of coordinates, or as the model zero based
            [(row, column)] location.
        modelgrid : flopy.discretization.StructuredGrid
            structured grid instance from flopy
        fmt : str
            format of the pour point information. Acceptable types include
            "xy", "rowcol", and "shp"

        Returns
        -------
        np.ndarray : delinated array of the watershed (1 is active,
        0 is inactive)
        """

        # get and store the watershed pour point for later cascade calcs.
        cwpp = self._get_pour_points(pour_point, modelgrid, fmt)[0]
        wpp = self._shape[1] + (self._shape[1] * cwpp[0]) + cwpp[1] + 1
        self._wpp = wpp

        # define the watershed
        watershed = self.define_subbasins(pour_point, modelgrid, fmt=fmt)

        # update internal flow direction and hrutype?
        watershed = np.pad(watershed, 1, "constant", constant_values=0).ravel()

        self._flow_directions[watershed == 0] = -1
        self._hru_type[watershed == 0] = 0

        watershed.shape = self._shape
        return watershed[1:-1, 1:-1]

    def _get_pour_points(self, pour_points, modelgrid, fmt):
        """
        Method to get the pour point locations

        Parameters
        ----------
        pour_points : str, tuple
            pour point can be supplied as either a shapefile, a list with an
            [(x, y)] tuple of coordinates, or as the model zero based
            [(row, column)] location.
        modelgrid : flopy.discretization.StructuredGrid
            structured grid instance from flopy
        fmt : str
            format of the pour point information. Acceptable types include
            "xy", "rowcol", and "shp"

        Returns
        -------
            list of row column tuples
        """
        if fmt == "rowcol":
            pour_points = pour_points
        elif fmt == "xy":
            pour_points = [
                modelgrid.intersect(pt[0], pt[1]) for pt in pour_points
            ]
        elif fmt == "shp":
            shape = shapefile.Reader(pour_points)
            pour_points = [pt.points[0] for pt in shape.shapes()]
            pour_points = [
                modelgrid.intersect(pt[0], pt[1]) for pt in pour_points
            ]
        else:
            raise Exception("pour point format not recognized.")
        return pour_points

    def define_subbasins(self, pour_points, modelgrid, fmt="rowcol"):
        """
        Method to perform sub-basin delineation on the flow direction array

        Parameters
        ----------
        pour_points : str, tuple
            pour point can be supplied as either a shapefile, a list with an
            [(x, y)] tuple of coordinates, or as the model zero based
            [(row, column)] location.
        modelgrid : flopy.discretization.StructuredGrid
            structured grid instance from flopy
        fmt : str
            format of the pour point information. Acceptable types include
            "xy", "rowcol", and "shp"

        Returns
        -------
        np.ndarray : delinated array of the watershed (1 is active,
        0 is inactive)
        """
        pour_points = self._get_pour_points(pour_points, modelgrid, fmt)

        # create pour point array
        pp_array = np.zeros((modelgrid.nrow, modelgrid.ncol), dtype=int)
        pp_list = []

        for index, (row, col) in enumerate(pour_points):
            ppnum = index + 1
            pp_array[(row, col)] = ppnum
            pp_list.append(ppnum)

        # identify flow directions
        if self._closed_basin:
            flow_directions = self._cb_flow_directions
        else:
            flow_directions = self._flow_directions

        # format pour point array to work with flow directions
        # pad and convert to 1d array
        pp_array = np.pad(pp_array, 1, "constant", constant_values=0).ravel()

        # get each pour points start index
        pp_idxs = {p: np.argwhere(pp_array == p)[0][0] for p in pp_list}

        # use flow direction to build pour point connection graph
        pp_conns = {p: {"from": [], "to": []} for p in pp_list}
        # loop through pour points
        for pp in pp_list:
            # get pout point index
            pp_idx = pp_idxs[pp]
            # pour point cell flow direction
            pp_fdir = flow_directions[pp_idx]
            # get the next cell
            n = pp_idx + self._offset_dict[pp_fdir]
            # loop through flow direction array till dead end
            # or another pour point
            while self._next_cell(n):
                # check flow direction
                ndir = flow_directions[n]
                if np.isnan(ndir):
                    break

                # check self intersection
                if pp_array[n] == pp:
                    break

                # check if there is a new pour point
                if pp_array[n] != 0 and pp_array[n] != pp:
                    # update pour point connection graph
                    pp_conns[pp_array[n]]["from"].append(pp)
                    pp_conns[pp]["to"].append(pp_array[n])
                    break

                # update pour point array path
                pp_array[n] = pp

                # get next cell
                n = n + self._offset_dict[ndir]

        # build order stack for pour points based on connection graph
        # identify "outlet" points
        outlets = [n for n in pp_conns if len(pp_conns[n]["to"]) == 0]
        # set up and order stack and processing stack
        order_stack = []
        proc_stack = []
        # loop through outlets
        for op in outlets:
            # append the outlet pour point to order stack and processing stack
            order_stack.append(op)
            proc_stack.append(op)
            # process rest of tree
            while len(proc_stack) > 0:
                proc = proc_stack[0]
                if proc not in order_stack:
                    order_stack.append(proc)
                proc_stack.extend(pp_conns[proc]["from"])
                proc_stack.remove(proc)

        order_stack = list(reversed(order_stack))
        # intialize watershed array as zeros
        subbasins = np.zeros(flow_directions.shape, dtype=int)
        # for ix, pp_id in pour_point_dict.items():
        for pp in order_stack:
            ix = pp_idxs[pp]
            subbasins[ix] = pp
            idxs = self._offsets + ix
            # get indexes of neighbor cells that flow into cell
            flow_idxs = idxs[
                (flow_directions[idxs] == self._d8_vectors_r)
            ].tolist()
            stack = flow_idxs
            while stack:
                for ixx in stack:
                    if subbasins[ixx] != 0:
                        stack = stack[1:]
                        continue
                    stack = stack[1:]
                    subbasins[ixx] = pp
                    idxs = self._offsets + ixx
                    flow_idxs = idxs[
                        (flow_directions[idxs] == self._d8_vectors_r)
                    ].tolist()
                    stack += flow_idxs

        subbasins.shape = self._shape
        return subbasins[1:-1, 1:-1]

    def _build_stream_nidp(self, streams_ix, fdir_array):
        """
        Method to build stream number of input drainage paths

        Parameters
        ----------
        streams_ix : list
            list of stream indexes
        fdir_array : np.ndarray
            flow direction array
        Returns
        -------
            np.ndarray
        """
        nidp_array = np.zeros(self._flow_directions.shape)
        nidp_dir = [2, 4, 8, 1, 16, 128, 64, 32]
        # loop through inners
        for ix in streams_ix:
            # get neighbors
            if np.isnan(fdir_array[ix]):
                raise Exception("Flow Direction not calculated")
            idxs = self._offsets + ix
            idxs_dir = fdir_array[idxs]
            nidp = sum(
                [
                    1
                    for nidpdir, idxdir in zip(nidp_dir, idxs_dir)
                    if nidpdir == idxdir
                ]
            )
            nidp_array[ix] = nidp
        return nidp_array

    def make_streams(
        self,
        fdir_array,
        fa_array,
        threshold,
        min_stream_len=None,
        max_reach=None,
        default_slope=0.001,
        min_slope=0.0001,
        max_slope=1.0,
    ):
        """
        Method to make streams for both PRMS and MODFLOW SFR

        Parameters
        ----------
        fdir_array : np.ndarray
            flow direction array
        fa_array : np.ndarray
            flow accumulation array
        threshold : float
            accumulated area threshold for defining streams
        min_stream_len : int
            optional minimum stream length in number of cells
        max_reach : int
            optional maximum number of reach cells per segment
        default_slope : float
            default slope value for setting SFR slope of all cells
        min_slope : float
            minimum slope calculated by make streams for SFR
        max_slope : float
            maximum SFR slope calculated by make streams

        Returns
        -------
            _StreamsObj : object containing SFR and PRMS stream information
        """
        # check and make flow direction, and accumulation
        assert fdir_array.shape == fa_array.shape

        # assign streams that are >= to threshold stream array has
        # flow accumulation stored
        streams = fa_array.copy()
        streams = np.where(streams > threshold, streams, np.nan)
        # get stream indexes
        streams = np.pad(
            streams, 1, "constant", constant_values=np.nan
        ).ravel()
        streams_ix = np.argwhere(~np.isnan(streams))

        # set up flow direction
        fdir_array = np.pad(
            fdir_array, 1, "constant", constant_values=np.nan
        ).ravel()
        streams_fdir_array = np.copy(fdir_array)

        # set fdir_array to nan where streams nan
        streams_fdir_array[np.isnan(streams)] = np.nan

        # build nidp for streams
        nidp = self._build_stream_nidp(streams_ix, streams_fdir_array)

        # find headwaters (starts of streams) indexes of possible headwaters
        options = np.argwhere(~np.isnan(streams))
        headwaters = np.zeros_like(fdir_array)

        # loop through possible headwaters
        for ix in options:
            if self._hru_type[ix] != 1:
                continue
            # get neighbors
            idxs = self._offsets + ix
            # get the fdirs
            inflow = streams_fdir_array[idxs] == self._d8_vectors_r
            if not inflow.any():
                headwaters[ix] = 1

        headwaters = np.argwhere(headwaters == 1).ravel()

        # segment array
        isegs = np.zeros(fdir_array.shape, dtype=int)

        # reach array ( array for now may change)
        reaches = np.zeros(fdir_array.shape, dtype=int)

        # outseg
        outsegs = np.zeros(fdir_array.shape, dtype=int)

        # dictionary with stream information
        stream_dict = {}
        iseg = 0
        for ix in headwaters:
            # set cell; n = r, c
            sstack = []
            n = ix
            iseg += 1
            reach = 0

            # trace down stream till there is no cell
            proc = True
            stream_dict[iseg] = {"inseg": [], "graph": []}
            while self._next_cell(n):
                if n in sstack:
                    break

                sstack.append(n)
                # check flow direction
                ndir = fdir_array[n]
                if np.isnan(ndir):
                    break

                # check if intersection cell
                if nidp[n] >= 2:
                    # check if cell has already been id
                    if isegs[n] != 0:
                        outseg = isegs[n]
                        outsegs[isegs == iseg] = outseg
                        stream_dict[iseg]["outseg"] = outseg
                        stream_dict[isegs[n]]["inseg"].append(iseg)
                        proc = False
                        break

                    # setting outseg
                    outsegs[isegs == iseg] = iseg + 1
                    stream_dict[iseg]["outseg"] = iseg + 1
                    iseg += 1

                    # set inseg
                    stream_dict[iseg] = {"inseg": [iseg - 1], "graph": []}

                    reach = 0

                # if reach = max reach reset reach and create new seg
                if reach == max_reach:
                    if iseg not in stream_dict:
                        raise Exception("iseg not in stream dict")
                    else:
                        if "outseg" in stream_dict[iseg]:
                            raise Exception("outseg already defined")
                        stream_dict[iseg]["outseg"] = iseg + 1
                        outsegs[isegs == iseg] = iseg + 1
                    reach = 0
                    iseg += 1
                    stream_dict[iseg] = {"inseg": [iseg - 1], "graph": []}

                reach += 1

                # set the segnum
                isegs[n] = iseg
                reaches[n] = reach
                stream_dict[iseg]["graph"].append(n)

                n = n + self._offset_dict[ndir]

            if proc:
                if iseg not in stream_dict:
                    raise Exception("iseg not in stream dict")
                else:
                    if "outseg" in stream_dict[iseg]:
                        raise Exception("outseg already defined")
                    stream_dict[iseg]["outseg"] = 0

        # remove streams less than min length
        unique, counts = np.unique(isegs, return_counts=True)
        # make list of pass through cells we do not want to remove
        dont_remove = []
        for seg in stream_dict:
            if (
                len(stream_dict[seg]["inseg"]) > 0
                and stream_dict[seg]["outseg"] != 0
            ):
                dont_remove.append(seg)

        if min_stream_len is not None:

            # get isegs to remove
            remove = [
                x
                for x, count in zip(unique, counts)
                if count < min_stream_len and x not in dont_remove
            ]

            # remove isegs from arrays
            for ix in remove:
                reaches[isegs == ix] = 0
                outsegs[isegs == ix] = 0
                outsegs[outsegs == ix] = 0
                isegs[isegs == ix] = 0
                del stream_dict[ix]

            # update outsegs in stream_dict as well
            for seg in stream_dict:

                if stream_dict[seg]["outseg"] in remove:
                    stream_dict[seg]["outseg"] = 0
                # update insegs
                stream_dict[seg]["inseg"] = [
                    x for x in stream_dict[seg]["inseg"] if x not in remove
                ]

        # renumber streams by topological order
        topology = Topology()
        for iseg, rec in stream_dict.items():
            if iseg > 0:
                try:
                    assert rec["outseg"] != -1
                except AssertionError:
                    raise AssertionError
                topology.add_connection(iseg, rec["outseg"])

        stack = topology.sort()

        renum_dict = {oldseg: ix + 1 for ix, oldseg in enumerate(stack)}
        renum_dict[0] = 0

        # renumber array
        isegs_new = np.zeros_like(isegs)
        outsegs_new = np.zeros_like(outsegs)
        for oldseg, newseg in sorted(renum_dict.items()):  # check this
            # renumber arrays
            isegs_new[isegs == oldseg] = newseg
            outsegs_new[outsegs == oldseg] = newseg

        isegs = isegs_new
        outsegs = outsegs_new

        # renumber stream dict
        new_stream_dict = {}
        for oldseg, stream_info in stream_dict.items():
            newseg = renum_dict[oldseg]
            new_stream_dict[newseg] = {
                "graph": stream_info["graph"],
                "outseg": renum_dict[stream_info["outseg"]],
                "inseg": [renum_dict[x] for x in stream_info["inseg"]],
            }

        stream_dict = new_stream_dict

        # build sfrtop from dem
        strtop = np.zeros(self._data.shape) * np.nan
        strtop[sstack] = self._data[sstack]

        # intialize rchlen and slope arrays
        rchlens = np.zeros(fdir_array.shape)
        slopes = np.zeros(fdir_array.shape)

        for seg in stream_dict:
            graph = np.array(stream_dict[seg]["graph"])
            ix = graph[0]
            insegs = stream_dict[seg]["inseg"]
            outseg = stream_dict[seg]["outseg"]

            xcenters = self._xcenters[graph]
            ycenters = self._ycenters[graph]

            xdiff1 = abs(xcenters[1:] - xcenters[:-1]) ** 2
            xdiff2 = abs(xcenters[:-1] - xcenters[1:]) ** 2

            ydiff1 = abs(ycenters[1:] - ycenters[:-1]) ** 2
            ydiff2 = abs(ycenters[:-1] - ycenters[1:]) ** 2

            xdiff1 = np.insert(xdiff1, 0, 0)
            ydiff1 = np.insert(ydiff1, 0, 0)

            xdiff2 = np.append(xdiff2, 0)
            ydiff2 = np.append(ydiff2, 0)

            dist1 = np.sqrt(xdiff1 + ydiff1) / 2.0
            dist2 = np.sqrt(xdiff2 + ydiff2) / 2.0

            dist = dist1 + dist2

            rchlens[graph] += dist

            # process insegs
            if len(insegs) > 0:
                inseg_ixs = np.array(
                    [stream_dict[inseg]["graph"][-1] for inseg in insegs]
                )

                xcenters = self._xcenters[inseg_ixs]
                ycenters = self._ycenters[inseg_ixs]

                xdiff = abs(xcenters - self._xcenters[ix]) ** 2
                ydiff = abs(ycenters - self._ycenters[ix]) ** 2

                dist = np.sqrt(xdiff + ydiff) / 2.0

                rchlens[ix] += dist.sum()

            # process outseg
            if outseg != 0:
                outseg_ix = stream_dict[outseg]["graph"][0]

                xcenters = self._xcenters[outseg_ix]
                ycenters = self._ycenters[outseg_ix]

                xdiff = abs(xcenters - self._xcenters[ix]) ** 2
                ydiff = abs(ycenters - self._ycenters[ix]) ** 2

                dist = np.sqrt(xdiff + ydiff) / 2.0
                rchlens[ix] += dist

            else:
                # get fdir, xcenters process
                fdir = fdir_array[ix]
                idx = ix + self._offset_dict[fdir]

                xdiff = abs(self._xcenters[idx] - self._xcenters[ix])
                ydiff = abs(self._ycenters[idx] - self._ycenters[ix])

                dist = np.sqrt(xdiff + ydiff) / 2.0
                rchlens[ix] += dist

        # calculate slope
        for seg in stream_dict:
            stream_cells = stream_dict[seg]["graph"]
            for idx, cell in enumerate(stream_cells):
                if (idx + 1) < len(stream_cells):
                    next_cell = stream_cells[idx + 1]
                    rchlen = rchlens[cell]
                    nslope = (strtop[cell] - strtop[next_cell]) / rchlen
                    slopes[cell] = nslope
                else:
                    slopes[cell] = default_slope

        slopes[slopes < min_slope] = min_slope
        slopes[slopes > max_slope] = max_slope

        unique, counts = np.unique(isegs, return_counts=True)
        max_reaches = {k: v for k, v in zip(unique, counts) if k != 0}

        # irunbnd
        irunbnd = np.zeros_like(isegs)
        for iseg in unique:
            if iseg != 0:
                proc = True
                ixs = np.argwhere(isegs == iseg)
                irunbnd[ixs] = iseg
                while proc:

                    idxs = self._offsets + ixs
                    idx_flow = fdir_array[idxs]
                    inflow = idx_flow == self._d8_vectors_r
                    ixxs = idxs[inflow]
                    iseg_check = isegs[ixxs]
                    irun_check = irunbnd[ixxs]
                    check = np.where(
                        (iseg_check == 0) & (irun_check == 0), True, False
                    )
                    ixxs = ixxs[check]
                    if len(ixxs) == 0:
                        proc = False
                        break
                    irunbnd[ixxs] = iseg
                    ixs = np.reshape(ixxs, (len(ixxs), 1))

        # set to 2d for output
        irunbnd.shape = self._shape
        irunbnd = irunbnd[1:-1, 1:-1]

        isegs.shape = self._shape
        isegs = isegs[1:-1, 1:-1]

        reaches.shape = self._shape
        reaches = reaches[1:-1, 1:-1]

        rchlens.shape = self._shape
        rchlens = rchlens[1:-1, 1:-1]

        slopes.shape = self._shape
        slopes = slopes[1:-1, 1:-1]

        strtop.shape = self._shape
        strtop = strtop[1:-1, 1:-1]

        # # build reach data
        reach_data_dict = {}
        for iseg in list(sorted(unique)):
            if iseg != 0:
                if iseg not in reach_data_dict:
                    reach_data_dict[iseg] = {}
                ixs = np.argwhere(isegs == iseg)
                for ix in ixs:
                    ireach = reaches[ix[0], ix[1]]
                    reach_data_dict[iseg][ireach] = {}
                    reach_data_dict[iseg][ireach]["i"] = ix[0]
                    reach_data_dict[iseg][ireach]["j"] = ix[1]
                    reach_data_dict[iseg][ireach]["rchlen"] = rchlens[
                        ix[0], ix[1]
                    ]
                    reach_data_dict[iseg][ireach]["strtop"] = strtop[
                        ix[0], ix[1]
                    ]
                    reach_data_dict[iseg][ireach]["slope"] = slopes[
                        ix[0], ix[1]
                    ]

        nreaches = np.count_nonzero(reaches)
        reach_data = ModflowSfr2.get_empty_reach_data(nreaches)
        reach_data_iseg = []
        reach_data_ireach = []
        for iseg, max_reach in sorted(max_reaches.items()):
            reach_data_iseg.extend([iseg] * max_reach)
            reach_data_ireach.extend([r for r in range(1, max_reach + 1)])

        reach_data["iseg"] = reach_data_iseg
        reach_data["ireach"] = reach_data_ireach
        reach_data["i"] = [
            reach_data_dict[seg][reach]["i"]
            for seg, reach in zip(reach_data_iseg, reach_data_ireach)
        ]
        reach_data["j"] = [
            reach_data_dict[seg][reach]["j"]
            for seg, reach in zip(reach_data_iseg, reach_data_ireach)
        ]
        reach_data["rchlen"] = [
            reach_data_dict[seg][reach]["rchlen"]
            for seg, reach in zip(reach_data_iseg, reach_data_ireach)
        ]
        reach_data["strtop"] = [
            reach_data_dict[seg][reach]["strtop"]
            for seg, reach in zip(reach_data_iseg, reach_data_ireach)
        ]
        reach_data["slope"] = [
            reach_data_dict[seg][reach]["slope"]
            for seg, reach in zip(reach_data_iseg, reach_data_ireach)
        ]

        # apply sfr reach defaults reach_defaults = bd.defaults['sfr']['reach']
        reach_defaults = self._defaults["sfr"]["reach"]
        for key, val in reach_defaults.items():
            val = [val] * nreaches
            reach_data[key] = val

        # build the segment_data array
        nsegments = len(stream_dict.keys())
        segment_data = ModflowSfr2.get_empty_segment_data(nsegments)
        segment_data["nseg"] = list(sorted(stream_dict.keys()))
        segment_data["outseg"] = [
            v["outseg"] for k, v in sorted(stream_dict.items())
        ]

        # apply segment defaults
        segment_defaults = self._defaults["sfr"]["segment"]
        for key, val in segment_defaults.items():
            val = [val] * nsegments
            segment_data[key] = val

        # make pour points
        pour_points = np.zeros(self._size)
        ppct = 0
        for seg in stream_dict:
            if stream_dict[seg]["outseg"] == 0:
                end_cell = stream_dict[seg]["graph"][-1]
                ppct += 1
                pour_points[end_cell] = 1

        pour_points.shape = self._shape
        pour_points = pour_points[1:-1, 1:-1]

        # calculate aspect (using the fdir??)
        aspect = np.zeros_like(fdir_array)
        aspect[fdir_array == 1] = 0
        aspect[fdir_array == 128] = 45
        aspect[fdir_array == 64] = 90
        aspect[fdir_array == 32] = 135
        aspect[fdir_array == 16] = 180
        aspect[fdir_array == 8] = 225
        aspect[fdir_array == 4] = 270
        aspect[fdir_array == 2] = 315

        aspect.shape = self._shape
        aspect = aspect[1:-1, 1:-1]

        gridded_data = {
            "iseg": isegs,
            "ireach": reaches,
            "outseg": outsegs,
            "irunbnd": irunbnd,
            "sfrtop": strtop,
            "rchlen": rchlens,
            "slope": slopes,
            "aspect": aspect,
            "pourpoints": pour_points,
        }

        streams_obj = _StreamsObj(reach_data, segment_data, gridded_data)
        return streams_obj

    @staticmethod
    def load_streams(f):
        """
        Method to load a saved binary streams file

        Parameters
        ----------
        f : str
            file name

        Returns
        -------
            _StreamsObj
        """
        return _StreamsObj.load(f)

    @staticmethod
    def load_cascades(f):
        """
        Method to load a saved binary cascades file

        Parameters
        ----------
        f : str
            file name

        Returns
        -------
            _Cascades
        """
        return _Cascades.load(f)


class _StreamsObj(object):
    """
    Class to hold information about the streams created with the
    FlowAccumulation tool

    Parameters
    ----------
    reach_data : np.recarray
        record array that corresponds to flopy SFR1 reach_data array
    segment_data : np.recarray
        record array the corresponds to flopy Sfr1 segment_data array
    gridded_data : dict
        dictionary of gridded data that can be later used for PRMS applications

    """

    def __init__(self, reach_data, segment_data, gridded_data):
        self.iseg = gridded_data["iseg"]
        self.ireach = gridded_data["ireach"]
        self.outseg = gridded_data["outseg"]
        self.irunbnd = gridded_data["irunbnd"]
        self.sfrtop = gridded_data["sfrtop"]
        self.rchlen = gridded_data["rchlen"]
        self.slope = gridded_data["slope"]
        self.aspect = gridded_data["aspect"]
        self.pour_points = gridded_data["pourpoints"]
        self.reach_data = reach_data
        self.segment_data = segment_data

    def write(self, f):
        """
        Method to write a binary streamsobj file

        Parameters
        ----------
        f : str
            file name

        Returns
        -------
            None
        """
        gsflow_io._write_pickle(f, self)

    @staticmethod
    def load(f):
        """
        Method to load a binary file containing a _StreamsObj object

        Parameters
        ----------
        f : str
            file name

        Returns
        -------
            _StreamsObj
        """
        return gsflow_io._read_pickle(f)


class _Cascades(object):
    """
    Object to hold Cascade results for prms

    Parameters
    ----------
    hru_up_id : np.ndarray
        array of hru_up_ids
    hru_down_id : np.ndarray
        array of hru_down_ids
    hru_pct_up : np.ndarray
        array of percentage of flow from up hru
    hru_strmseg_down_id : np.ndarray
        array of stream seg id's a cascade connects to
    """

    def __init__(
        self, hru_up_id, hru_down_id, hru_pct_up, hru_strmseg_down_id=None
    ):
        self.ncascade = hru_up_id.size
        self.hru_up_id = hru_up_id
        self.hru_down_id = hru_down_id
        self.hru_pct_up = hru_pct_up
        self.hru_strmseg_down_id = hru_strmseg_down_id

    def write(self, f):
        """
        Method to write a binary cascades file

        Parameters
        ----------
        f : str
            file name

        Returns
        -------
            None
        """
        gsflow_io._write_pickle(f, self)

    @staticmethod
    def load(f):
        """
        Method to load a binary file containing a _Cascades object

        Parameters
        ----------
        f : str
            file name

        Returns
        -------
            _Cascades
        """
        return gsflow_io._read_pickle(f)

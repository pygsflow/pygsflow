import numpy as np
import os
import warnings
import flopy as fp

try:
    import shapefile
except ImportError:
    shapefile = None

try:
    import pycrs
except ImportError:
    pycrs = None

warnings.simplefilter("always", UserWarning)


class Modsim(object):
    """
    Class to handle creating MODSIM inputs
    for the MODSIM-GSFLOW coupled modeling option
    in GSFLOW.

    Parameters
    ----------
        model : gsflow.GsflowModel instance or gsflow.modflow.Modflow instance
        other : str or None
            include can be used to add an optional sfr field (ex. strhc1) to the
            stream vector shapefile.

    Examples
    --------

    >>> gsf = gsflow.GsflowModel.load_from_file("gsflow.control")
    >>> modsim = gsflow.modsim.Modsim(gsf)
    >>> modsim.write_modsim_shapefile("myshp.shp")

    """

    def __init__(self, model, other=None):
        from ..gsflow import GsflowModel
        from ..modflow import Modflow

        if isinstance(model, (Modflow, fp.modflow.Modflow)):
            self.parent = model
            self.mf = model
        else:
            self.parent = model
            self.mf = self.parent.mf
        if other is not None:
            other = other.lower()
        self._other = other
        self._sfr = self.mf.get_package("SFR")
        self._lak = self.mf.get_package("LAK")
        self._ready = True

        if self.mf is None:
            self._ready = False

        if self._sfr is None:
            self._ready = False

        self._nearest = True

    @property
    def sfr_segs(self):
        """
        Get all of the SFR segments

        Returns
        -------
            list of sfr segments

        """
        if not self._ready:
            return ()
        else:
            reach_data = self._sfr.reach_data
            segs = []
            for rec in reach_data:
                segs.append(rec.iseg)

            segs = list(set(segs))
            return sorted(segs)

    @property
    def lake_segs(self):
        """
        Get all of the lake numbers connected to SFR segments

        Returns
        -------
            list of all lakes in the sfr network

        """
        if not self._ready:
            return ()
        elif self._lak is None:
            return ()
        else:
            seg_data = self._sfr.segment_data
            lakes = []
            for per, recarray in seg_data.items():
                for rec in recarray:
                    if rec.iupseg < 0:
                        lakes.append(rec.iupseg)
                    if rec.outseg < 0:
                        lakes.append(rec.outseg)

            lakes = list(set(lakes))
            return sorted(lakes)

    @property
    def sfr_topology(self):
        """
        Creates a list of SFR topology objects

        Returns
        -------
            list of _SfrTopology objects for writing
            to shapefile

        """
        sfr_topo = []
        for seg in self.sfr_segs:
            sfr_topo.append(_SfrTopology(self._sfr, self.mf, seg, self._other))

        return sfr_topo

    @property
    def lake_topology(self):
        """
        Creates a list of LAK topology objects

        Returns
        -------
            list of _LakTopology objects for writing
            to shapefile

        """
        lake_topo = []
        for lake in self.lake_segs:
            lake_topo.append(
                _LakTopology(self._lak, self.mf, lake, self._nearest)
            )

        return lake_topo

    def __set_spillway_flag(self, flag, sfr_topology):
        """
        Flips the sfr topology flag based on user supplied
        flag.

        Parameters
        ----------
        flag : str or list
            "elev" for elevation rules
            "flow" for flow rules
            list of segments to flip certain segments
        sfr_topology : list
            list of _SfrTopology instances

        Returns
        -------
        sfr_topology : list
            list of _SfrTopology instances

        """
        temp = []
        if isinstance(flag, list):
            for sfr in sfr_topology:
                if sfr.attributes.iseg in flag:
                    sfr.attributes.spill_flg = 1
                else:
                    pass
                temp.append(sfr)

        else:
            res_connect = []
            for ix, sfr in enumerate(sfr_topology):
                if sfr.attributes.iupseg < 0:
                    res_connect.append(ix)

            t = []
            for i in res_connect:
                sfr = sfr_topology[i]

                for j in res_connect:
                    sfr2 = sfr_topology[j]
                    if sfr.attributes.iseg == sfr2.attributes.iseg:
                        continue
                    elif sfr.attributes.iupseg == sfr2.attributes.iupseg:
                        if sfr.attributes.outseg == sfr2.attributes.outseg:
                            if flag.lower() == "elev":
                                if sfr.attributes.elev > sfr2.attributes.elev:
                                    sfr.attributes.spill_flg = 1
                                else:
                                    pass
                            if flag.lower() == "flow":
                                if sfr.attributes.flow == 0:
                                    sfr.attributes.spill_flg = 1
                                else:
                                    pass

                t.append(sfr)

            for ix, sfr in enumerate(sfr_topology):
                if ix in res_connect:
                    temp.append(t.pop(0))
                else:
                    temp.append(sfr)

        return temp

    def __add_field(self, field_name):
        """

        Parameters
        ----------
        field_name

        Returns
        -------

        """

    def write_modsim_shapefile(
        self, shp=None, proj4=None, flag_spillway=False, nearest=True
    ):
        """
        Method to create a modsim compatible
        shapefile from GSFLOW model inputs (SFR, LAK)
        package.

        Parameters
        ----------
        shp : str
            optional shapefile name, if none
            will be written in gsflow directory using
            the model name.
        proj4 : str
            proj4 projection string, if none will try to
            grab proj4 or epsg from flopy modelgrid
        flag_spillway : bool, str, list
            if flag_spillway is indicated then MODSIM will change
            the spill_flg attribute to one. This can be accomplished
            by one of three methods.

            1.) flag_spillway="elev", the code will search for spillways
            from reservoirs based on elevation rules
            2.) flag_spillway="flow", the code will search for spillways
            from reservoirs based on flow rules
            3.) flag_spillway=[3, 4, 5, ...] a user supplied list of SFR
            segments can be supplied to flag spillways
        nearest : bool
            if nearest is True, lak topology will connect to the nearest
            SFR reaches based on the segment they are connected to. If False
            lak topology will connect to the start or end of the Segment based
            on iupseg and outseg

        """
        if not self._ready:
            return

        if shapefile is None:
            warn = "Pyshp must be installed to write MODSIM shapefile"
            warnings.warn(warn)
            return

        if shp is None:
            ws = self.parent.control.model_dir
            t = self._sfr.file_name[0].split(".")
            t = ".".join(t[:-1])
            name = t + "_modsim.shp"
            shp = os.path.join(ws, name)

        self._nearest = nearest

        sfr_topology = self.sfr_topology
        lake_topology = self.lake_topology

        if flag_spillway:
            sfr_topology = self.__set_spillway_flag(
                flag_spillway, sfr_topology
            )

        w = shapefile.Writer(shp)
        w.shapeType = 3
        w.field("ISEG", "N")
        w.field("IUPSEG", "N")
        w.field("OUTSEG", "N")
        w.field("SPILL_FLG", "N")
        if self._other is not None:
            w.field(self._other.upper(), "N", decimal=5)

        for sfr in sfr_topology:
            w.line(sfr.polyline)
            attributes = sfr.attributes
            if self._other is None:
                w.record(
                    attributes.iseg,
                    attributes.iupseg,
                    attributes.outseg,
                    attributes.spill_flg,
                )
            else:
                w.record(
                    attributes.iseg,
                    attributes.iupseg,
                    attributes.outseg,
                    attributes.spill_flg,
                    attributes.other,
                )

        for lake in lake_topology:
            for ix, attributes in enumerate(lake.attributes):
                w.line([lake.polyline[ix]])
                if self._other is None:
                    w.record(
                        attributes.iseg,
                        attributes.iupseg,
                        attributes.outseg,
                        attributes.spill_flg,
                    )
                else:
                    w.record(
                        attributes.iseg,
                        attributes.iupseg,
                        attributes.outseg,
                        attributes.spill_flg,
                        attributes.other,
                    )
        try:
            w.close()
        except AttributeError:
            pass

        if pycrs is None:
            warn = (
                "PyCRS must be installed to add a projection"
                " to {}".format(shp)
            )
            warnings.warn(warn)
            return

        t = shp.split(".")
        if len(t) == 1:
            prj = shp + ".prj"
        else:
            s = ".".join(t[:-1])
            prj = s + ".prj"

        if proj4 is None:
            proj4 = self.mf.modelgrid.proj4

        epsg = self.mf.modelgrid.epsg

        try:
            crs = pycrs.parse.from_proj4(proj4)
        except:
            crs = None

        if crs is None:
            try:
                crs = pycrs.parse.from_epsg_code(epsg)
            except:
                crs = None

        if crs is None:
            warn = (
                "Please provide a valid proj4 or epsg code to "
                "flopy's model grid: Skipping writing {}".format(prj)
            )
            warnings.warn(warn)
            return

        with open(prj, "w") as foo:
            foo.write(crs.to_esri_wkt())


class _LakTopology(object):
    """
    Object that creates lake centroids and
    defines the topology/connectivity of
    the LAK file to the SFR object

    Parameters
    ----------
    lak : flopy.modflow.ModflowLak object
    model : flopy.modflow.Modflow or gsflow.modflow.Modflow
    lakeno : int
        lake number
    nearest : bool
        defaults to True, creates a connection to nearest node
        if False creates connection to the last reach
    """

    def __init__(self, lak, model, lakeno, nearest=True):
        self._parent = model
        self._lak = lak
        self._sfr = self._parent.get_package("SFR")
        self._mg = self._parent.modelgrid
        self._xv = self._mg.xvertices
        self._yv = self._mg.yvertices
        self._lakeno = lakeno
        self._centroid = None
        self._connections = None
        self._polyline = None
        self._attributes = None
        self._ij = None
        self._nearest = nearest

    @property
    def lakeno(self):
        return self._lakeno

    @property
    def centroid(self):
        if self._centroid is None:
            self._set_lake_centroid()
        return self._centroid

    @property
    def connections(self):
        if self._connections is None:
            self._set_lake_connectivity()
        return self._connections

    @property
    def polyline(self):
        if self._polyline is None:
            self._set_polyline()
        return self._polyline

    @property
    def attributes(self):
        if self._attributes is None:
            self._set_lake_connectivity()
        return self._attributes

    def _set_lake_centroid(self):
        """
        Method to calculate the geometric
        centroid of a lake

        """
        # get a 3d array of lak locations
        lakes = self._lak.lakarr.array

        lakeno = self.lakeno
        if lakeno < 0:
            lakeno *= -1

        verts = []
        for lake3d in lakes:
            if len(lake3d.shape) == 3:
                t = np.where(lake3d == lakeno)
                i = list(t[1])
                j = list(t[2])
                ij = list(zip(i, j))
                for i, j in ij:
                    verts += self._mg.get_cell_vertices(i, j)

        if verts:
            verts = np.array(list(set(verts))).T
            xc = np.mean(verts[0])
            yc = np.mean(verts[1])
            self._centroid = (xc, yc)
        else:
            self._centroid = None

    def _set_lake_connectivity(self):
        """
        Method to define lake connections to
        the sfr network. Will be used to
        define polylines and topology

        """
        if self._sfr is None:
            return

        lakeno = self.lakeno
        if lakeno > 0:
            lakeno *= -1

        cseg = []
        attrs = []
        for per, recarray in self._sfr.segment_data.items():
            for rec in recarray:
                if rec.nseg in cseg:
                    continue
                else:
                    if rec.iupseg == lakeno:
                        cseg.append(rec.nseg)
                        attrs.append(_Attributes(lakeno, outseg=rec.nseg))
                    elif rec.outseg == lakeno:
                        cseg.append(rec.nseg)
                        attrs.append(_Attributes(lakeno, iupseg=rec.nseg))

        temp = []
        for seg in cseg:
            if seg not in temp:
                temp.append(seg)

        cseg = temp

        ij = []
        reach_data = self._sfr.reach_data
        reach_data.sort(axis=0, order=["iseg", "ireach"])
        for seg in cseg:
            temp = []
            for rec in reach_data:
                if rec.iseg == seg:
                    temp.append((rec.i, rec.j))

            ij.append(temp)

        # now we use the distance equation to connect to
        # the closest....
        if self._nearest:
            verts = []
            for conn in ij:
                tverts = []
                dist = []
                for i, j in conn:
                    xv = self._mg.xcellcenters[i, j]
                    yv = self._mg.ycellcenters[i, j]
                    tverts.append((xv, yv))
                    a = (xv - self.centroid[0]) ** 2
                    b = (yv - self.centroid[1]) ** 2
                    c = np.sqrt(a + b)
                    dist.append(c)

                if tverts:
                    vidx = dist.index(np.min(dist))
                    verts.append(tverts[vidx])
        else:
            verts = []
            for ix, conn in enumerate(ij):
                if attrs[ix].outseg == 0:
                    i, j = conn[-1]
                else:
                    i, j = conn[0]
                xv = self._mg.xcellcenters[i, j]
                yv = self._mg.ycellcenters[i, j]
                verts.append((xv, yv))

        if verts:
            self._connections = verts
            self._attributes = attrs

        else:
            self._connections = None
            self._attributes = None

    def _set_polyline(self):
        """
        Method to create and set the polyline input
        for pyshp

        """
        if self.connections is None:
            return
        if self.centroid is None:
            return

        polyline = []
        for ix, conn in enumerate(self.connections):
            if self.attributes[ix].iupseg == 0:
                part = [list(self.centroid), list(conn)]
            else:
                part = [list(conn), list(self.centroid)]

            polyline.append(part)

        self._polyline = polyline


class _SfrTopology(object):
    """
    Object that creates lake centroids and
    defines the topology/connectivity of
    a sfr segment

    Parameters
    ----------
    sfr : flopy.modflow.ModflowSfr object
    model : flopy.modflow.Modflow or gsflow.modflow.Modflow object
    iseg : int
        sfr segment number
    other : str or None
        include can be used to add an optional sfr field (ex. strhc1) to the
        stream vector shapefile.


    """

    def __init__(self, sfr, model, iseg, other=None):
        self._sfr = sfr
        self._parent = model
        self._iseg = iseg
        self._other = other
        self._mg = self._parent.modelgrid
        self._xv = self._mg.xvertices
        self._yv = self._mg.yvertices
        self._connections = None
        self._polyline = None
        self._attributes = None
        self._ij = None

    @property
    def iseg(self):
        return self._iseg

    @property
    def ij(self):
        if self._ij is None:
            self._set_attributes()
        return self._ij

    @property
    def connections(self):
        if self._connections is None:
            self._set_sfr_connectivity()
        return self._connections

    @property
    def polyline(self):
        if self._polyline is None:
            self._set_polyline()
        return [self._polyline]

    @property
    def attributes(self):
        if self._attributes is None:
            self._set_attributes()
        return self._attributes

    def _set_sfr_connectivity(self):
        """
        Method to get and set the the sfr
        connectivity and get the centroid of
        the upstream or downstream segment
        for shapefile exporting

        """
        if self._sfr is None:
            return

        outseg = self.attributes.outseg
        iupseg = self.attributes.iupseg

        ijup = []
        ijout = []
        for rec in self._sfr.reach_data:
            if rec.iseg == iupseg:
                ijup.append([rec.i, rec.j])
            elif rec.iseg == outseg:
                ijout.append([rec.i, rec.j])
            else:
                pass

        updist = []
        for i, j in ijup:
            a = (i - self.ij[0]) ** 2
            b = (j - self.ij[1]) ** 2
            c = np.sqrt(a + b)
            updist.append(c)

        outdist = []
        for i, j in ijout:
            a = (i - self.ij[0]) ** 2
            b = (j - self.ij[1]) ** 2
            c = np.sqrt(a + b)
            outdist.append(c)

        if updist:
            upidx = updist.index(np.min(updist))
            ijup = [ijup[upidx]]

        if outdist:
            outix = outdist.index(np.min(outdist))
            ijout = [ijout[outix]]

        vup = ()
        for i, j in ijup:
            xv = self._mg.xcellcenters[i, j]
            yv = self._mg.ycellcenters[i, j]
            vup = (xv, yv)

        vout = ()
        for i, j in ijout:
            xv = self._mg.xcellcenters[i, j]
            yv = self._mg.ycellcenters[i, j]
            vout = (xv, yv)

        self._connections = [vup, vout]

    def _set_polyline(self):
        """
        Method to construct a polyline for the
        SFR segments

        """
        if self._sfr is None:
            return

        conn = self.connections
        ij = []
        reach_data = self._sfr.reach_data
        reach_data.sort(axis=0, order=["iseg", "ireach"])
        for rec in reach_data:
            if rec.iseg == self.iseg:
                ij.append([rec.i, rec.j])

        line = []
        if conn[0]:
            line.append(list(conn[0]))

        for i, j in ij:
            xv = self._mg.xcellcenters[i, j]
            yv = self._mg.ycellcenters[i, j]
            line.append([xv, yv])

        if conn[-1]:
            line.append(list(conn[-1]))

        if line:
            self._polyline = line
        else:
            self._polyline = None

    def _set_attributes(self):
        """
        Method to set the attribute field
        for a SFR segment

        """
        if self._sfr is None:
            return

        outseg = []
        iupseg = []
        flow = []
        for per, recarray in self._sfr.segment_data.items():
            for rec in recarray:
                if rec.nseg == self.iseg:
                    outseg.append(rec.outseg)
                    iupseg.append(rec.iupseg)
                    flow.append(rec.flow)
                    break

        ij = []
        strtop = []
        other = []
        reach_data = self._sfr.reach_data
        reach_data.sort(axis=0, order=["iseg", "ireach"])
        for rec in reach_data:
            if rec.iseg == self.iseg:
                ij.append([rec.i, rec.j])
                strtop.append(rec.strtop)
                if self._other is not None:
                    other.append(rec[self._other])

        if ij:
            self._ij = tuple(ij[-1])

        if outseg:
            outseg = list(set(outseg))[0]
        else:
            outseg = 0

        if iupseg:
            iupseg = list(set(iupseg))[0]
        else:
            iupseg = 0

        strtop = min(strtop)
        if self._other is not None:
            other = np.mean(other)
        else:
            other = None
        flow = max(flow)

        self._attributes = _Attributes(
            self.iseg, iupseg, outseg, flow, strtop, other
        )


class _Attributes(object):
    """
    Object oriented storage method to standardize
    the topology function calls

    Parameters
    ----------
    iseg : int
        segment number
    iupseg : int
        upstream segment number
    outseg : int
        output segment number
    flow : float
        maximum specified flow rate
    strtop : float
        minimum strtop elevation
    other : None or float
        additional field data
    """

    def __init__(self, iseg, iupseg=0, outseg=0, flow=0, strtop=0, other=None):
        self.iseg = iseg
        self.iupseg = iupseg
        self.outseg = outseg
        self.flow = flow
        self.elev = strtop
        self.spill_flg = 0
        if other is None:
            self.other = np.nan
        else:
            self.other = other

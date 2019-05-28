import numpy as np
import os
import warnings

try:
    import shapefile
except ImportError:
    shapefile = None

try:
    import pycrs
except ImportError:
    pycrs = None

warnings.simplefilter('always', UserWarning)


class Modsim(object):
    """
    Class to handle creating MODSIM inputs
    for the MODSIM-GSFLOW coupled modeling option
    in GSFLOW.

    Parameters:
    ----------
        model : gsflow.GsFlowModel instance

    """
    def __init__(self, model):
        self.parent = model
        self.mf = self.parent.mf
        self._sfr = self.mf.get_package("SFR")
        self._lak = self.mf.get_package("LAK")
        self._ready = True

        if self.mf is None:
            self._ready = False

        if self._sfr is None:
            self._ready = False

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
            sfr_topo.append(_SfrTopology(self._sfr, self.mf, seg))

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
            lake_topo.append(_LakTopology(self._lak, self.mf, lake))

        return lake_topo

    def write_modsim_shapefile(self, shp=None, proj4=None):
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

        w = shapefile.Writer(shp)
        w.shapeType = 3
        w.field("ISEG", "N")
        w.field("IUPSEG", "N")
        w.field("OUTSEG", "N")

        for sfr in self.sfr_topology:
            w.line(sfr.polyline)
            attributes = sfr.attributes
            w.record(attributes['iseg'],
                     attributes['iupseg'],
                     attributes['outseg'])

        for lake in self.lake_topology:
            w.line(lake.polyline)
            attributes = lake.attributes
            w.record(attributes['iseg'],
                     attributes['iupseg'],
                     attributes['outseg'])
        try:
            w.close()
        except AttributeError:
            pass

        if pycrs is None:
            warn = "PyCRS must be installed to add a projection" \
                   " to {}".format(shp)
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
            warn = "Please provide a valid proj4 or epsg code to " \
                   "flopy's model grid: Skipping writing {}".format(prj)
            warnings.warn(warn)
            return

        with open(prj, "w") as foo:
            foo.write(crs.to_esri_wkt())


class _LakTopology(object):
    """
    Object that creates lake centroids and
    defines the topology/connectivity of
    the LAK file to the SFR object

    lak : flopy.modflow.ModflowLak object
    model : flopy.modflow.Modflow or gsflow.modflow.Modflow
    lakeno : int
        lake number

    """
    def __init__(self, lak, model, lakeno):
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
        self._ij = None

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
        lakeno = self.lakeno
        if lakeno > 0:
            lakeno *= -1
        return {"iseg": lakeno, "iupseg": 0, "outseg": 0}

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
        for per, recarray in self._sfr.segment_data.items():
            for rec in recarray:
                if rec.iupseg == lakeno:
                    cseg.append(rec.nseg)
                elif rec.outseg == lakeno:
                    cseg.append(rec.nseg)

        cseg = list(set(cseg))

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
        verts = []
        for conn in ij:
            tverts = []
            dist = []
            for i, j in conn:
                xv = self._mg.xcellcenters[i, j]
                yv = self._mg.ycellcenters[i, j]
                tverts.append((xv, yv))
                a = (xv - self.centroid[0])**2
                b = (yv - self.centroid[1])**2
                c = np.sqrt(a + b)
                dist.append(c)

            if tverts:
                vidx = dist.index(np.min(dist))
                verts.append(tverts[vidx])

        if verts:
            self._connections = verts

        else:
            self._connections = None

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
        for conn in self.connections:
            part = [list(self.centroid), list(conn)]
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

    """
    def __init__(self, sfr, model, iseg):
        self._sfr = sfr
        self._parent = model
        self._iseg = iseg
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

        outseg = self.attributes['outseg']
        iupseg = self.attributes['iupseg']

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
            a = (i - self.ij[0])**2
            b = (j - self.ij[1])**2
            c = np.sqrt(a + b)
            updist.append(c)

        outdist = []
        for i, j in ijout:
            a = (i - self.ij[0])**2
            b = (j - self.ij[1])**2
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
        for per, recarray in self._sfr.segment_data.items():
            for rec in recarray:
                if rec.nseg == self.iseg:
                    outseg.append(rec.outseg)
                    iupseg.append(rec.iupseg)
                    break

        ij = []
        reach_data = self._sfr.reach_data
        reach_data.sort(axis=0, order=["iseg", "ireach"])
        for rec in reach_data:
            if rec.iseg == self.iseg:
                ij.append([rec.i, rec.j])

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

        self._attributes = {"iseg": self.iseg,
                            "iupseg": iupseg ,
                            "outseg": outseg}

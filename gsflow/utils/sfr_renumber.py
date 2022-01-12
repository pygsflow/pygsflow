import numpy as np
import pandas as pd
import flopy as fp
import warnings

warnings.simplefilter("always", UserWarning)


class SfrRenumber(object):
    """
    Class method that facilitates SFR package renumbering based on
    topology trees, discretization elevations, STRTOP elevation, or a user defined
    scheme. This class can be used to renumber SFR file or it can
    be used to renumber all SFR associated files in a model ex. SFR, GAGE, AG...

    Parameters
    ----------
    model : <flopy.modflow.Modflow> object or <gsflow.modflow.Modflow> object
        Modflow model object from FloPy or pyGSFLOW
    sfr : <flopy.modflow.ModflowSfr1> object
        optional and not needed if model is provided
    dis : <flopy.modflow.ModflowDis> object
        optional and not needed if model is provided
    gage : <flopy.modflow.ModflowGage> object
        optional and not needed if model is provided
    lak : <flopy.modflow.ModflowLak> object
        optional and not needed if model is provided
    ag : <gsflow.modflow.ModflowAwu> object
        optional and not needed if <gsflow.modflow.Modflow>
        model is provided
    scheme : str
        renumbering scheme, can be "dis" to renumber by land
        surface elevation, "sfr" to renumber by stream bed elevation
        or "user" to renumber by a user supplied scheme.
    user_scheme : dict
        used only if "user" is specified for scheme. Dictionary
        of {current segment : new segment} numbers
        example
            {1: 6,
             2: 1,
             3: 2 ....}

    Examples
    --------
    Renumber a model based on the topology of the SFR network connections
    (Recommended method)

    >>> import flopy
    >>> ml = flopy.modflow.Modflow.load("sagehen.nam")
    >>> sfr_renumber = SfrRenumber(model=ml)  # default scheme="topology"
    >>> sfr_renumber.renumber_all()

    Renumber a model based on the elevation of TOP in the dis file

    >>> import flopy
    >>> ml = flopy.modflow.Modflow.load("sagehen.nam")
    >>> sfr_renumber = SfrRenumber(model=ml, scheme="dis")
    >>> sfr_renumber.renumber_all()

    Renumber a model based on the elevation of STRTOP in the sfr file

    >>> import gsflow
    >>> gsf = gsflow.GsflowModel.load_from_file('gsflow.control')
    >>> ml = gsf.ml
    >>> sfr_renumber = SfrRenumber(model=ml, scheme="sfr")
    >>> sfr_renumber.renumber_all()

    or

    >>> import flopy
    >>> ml = flopy.modflow.Modflow.load("sagehen.nam")
    >>> sfr_renumber = SfrRenumber(model=ml, scheme="sfr")
    >>> sfr_renumber.renumber_all()

    Renumber a model using a user defined renumbering scheme

    >>> import flopy
    >>> ml = flopy.modflow.Modflow.load("sagehen.nam")
    >>> # create and apply a segment inversion scheme
    >>> t = range(50, 0, -1)
    >>> my_scheme = {ix: i for ix, i in enumerate(t)}
    >>> sfr_renumber = SfrRenumber(model=ml, scheme="user", user_scheme=my_scheme)
    >>> # to renumber only SFR package
    >>> sfr_renumber.renumber_sfr()
    >>> # to renumber all packages that reference SFR segments
    >>> sfr_renumber.renumber_all()

    """

    def __init__(
        self,
        model=None,
        sfr=None,
        dis=None,
        gage=None,
        lak=None,
        ag=None,
        scheme="topology",
        user_scheme=None,
    ):
        from ..modflow import Modflow
        from ..gsflow import GsflowModel

        # check to see if user provided a GSflowModel instance
        if isinstance(model, GsflowModel):
            model = model.mf

        self.model = model
        if model is None:
            if sfr is None:
                err = "Either a model or SFR must be provided"
                raise AssertionError(err)

            if dis is None and scheme.lower() == "dis":
                err = "A Discretization object must be provided"
                raise AssertionError(err)

        self.sfr = sfr
        self.dis = dis
        self.gage = gage
        self.lak = lak
        self.ag = ag
        if model is not None:
            self.sfr = model.get_package("SFR")
            self.dis = model.get_package("DIS")
            self.gage = model.get_package("GAGE")
            self.lak = model.get_package("LAK")

            # check to see if the user provided a gsflow.modflow.Modflow object
            if isinstance(model, Modflow):
                self.ag = model.get_package("AWU")

        if self.sfr is None:
            raise AssertionError("SFR package not found, check model inputs")

        self._all_sfr_segments = None

        self.__scheme = scheme.lower()
        self.__renumber = None
        if scheme.lower() == "user":
            if user_scheme is None:
                err = "A numbering scheme must be provided"
                raise AssertionError(err)

            else:
                if 0 not in user_scheme:
                    user_scheme[0] = 0

                self.__renumber = user_scheme

    @property
    def all_sfr_segments(self):
        """
        Method to get a recarray of all SFR segments within a simulation,
        since all segments do not have to be active in a given stress
        period

        Returns
        -------
        ra : (np.recarray)
            recarray contains a single entry of segment data for each stream
            segment

        """
        if self._all_sfr_segments is None:
            ra = self.sfr.get_empty_segment_data(self.sfr.nss)
            i = 0
            for _, recarray in sorted(self.sfr.segment_data.items()):
                if i == self.sfr.nss:
                    break
                else:
                    for rec in recarray:
                        if rec.nseg in ra.nseg:
                            pass
                        else:
                            ra[i] = rec
                            i += 1

                        if i == self.sfr.nss:
                            break
            self._all_sfr_segments = ra

        return self._all_sfr_segments

    @property
    def renumbering(self):
        if self.__renumber is None:
            self.__calculate_renumbering_scheme()
        return self.__renumber

    @property
    def scheme(self):
        return self.__scheme

    def __calculate_renumbering_scheme(self):
        """
        Method to create the renumbering dictionary
        """
        nss = self.sfr.nss

        if self.scheme in ("sfr", "dis"):
            # sort by topography
            data = pd.DataFrame(columns=["segment", "elev"])
            sfr_ds2 = self.sfr.reach_data

            if self.scheme == "sfr":
                for iseg in range(1, nss + 1):
                    for record in sfr_ds2:
                        if record.iseg == iseg:
                            if record.ireach == 1:
                                strtop = record.strtop
                                s2 = pd.Series(
                                    {"segment": iseg, "elev": strtop}
                                )
                                data = data.append(s2, ignore_index=True)
                                break

            elif self.scheme == "dis":
                top = self.dis.top.array
                for iseg in range(1, nss + 1):
                    for record in sfr_ds2:
                        if record.iseg == iseg:
                            if record.ireach == 1:
                                i = record.i
                                j = record.j
                                strelev = top[i, j]
                                s2 = pd.Series(
                                    {"segment": iseg, "elev": strelev}
                                )
                                data = data.append(s2, ignore_index=True)
                                break

            data = data.sort_values(by=["elev"], ascending=False)

        else:
            # sort by topology, preferred method
            segments = self.all_sfr_segments
            topo = Topology(nss)

            chk = []
            for rec in segments:
                topo.add_connection(rec.nseg, rec.outseg)
                chk.append(rec.outseg)
                if rec.iupseg < 0:
                    topo.add_connection(rec.iupseg, rec.nseg)

            # check for dead end lake diversions
            for seg in chk:
                if seg not in topo.topology:
                    topo.add_connection(seg, 0)

            stack = topo.sort()
            stack = [[iseg] for iseg in stack if iseg != 0]
            data = pd.DataFrame(
                stack,
                columns=[
                    "segment",
                ],
            )

        if self.__renumber is None:
            # add zeros for ioutseg renumbering
            self.__renumber = {0: 0}

        iseg = 1
        for ix, record in data.iterrows():
            if record.values[0] < 0:
                # trap for lak numbers
                self.__renumber[int(record.values[0])] = int(record.values[0])
            else:
                self.__renumber[int(record.values[0])] = iseg
                iseg += 1

    def renumber_sfr(self):
        """
        Method to renumber only the Sfr Package. Using
        renumber_all() is recommended over this method
        in most cases
        """
        if self.gage is not None:
            warn = (
                "Gage package will not be renumbered\n"
                "MODFLOW may not run properly"
            )
            warnings.warn(warn, UserWarning)
        elif self.ag is not None:
            warn = (
                "Ag package will not be renumbered\n"
                "GSFLOW may not run properly"
            )
            warnings.warn(warn, UserWarning)

        self.__renumber_sfr()

    def renumber_all(self):
        """
        User method to renumber SFR and all SFR related packages
        attached to a modflow model
        """
        self.__renumber_sfr()
        self.__renumber_gage()
        self.__renumber_ag()

    def __renumber_sfr(self):
        """
        Resuable renumbering scheme for SFR, called
        by both renumber_sfr and renumber all
        """
        # first renumber dataset 2
        nstrm = abs(self.sfr.nstrm)
        sfr_ds2 = self.sfr.reach_data
        ra = np.recarray((nstrm,), dtype=sfr_ds2.dtype)

        for ix, record in enumerate(sfr_ds2):
            record.iseg = self.renumbering[record.iseg]
            ra[ix] = record

        ra.sort(order=["iseg", "ireach"])

        # set reach_data and set stress_period_data
        self.sfr.reach_data = ra
        self.sfr.stress_period_data = fp.utils.MfList(
            self.sfr, ra, dtype=ra.dtype
        )

        # then renumber transient data
        sfr_ds6 = self.sfr.segment_data
        spd = {}

        for kper, mflist in sfr_ds6.items():
            ntrans = len(mflist)
            ra = np.recarray((ntrans,), dtype=mflist.dtype)

            for ix, record in enumerate(mflist):
                # need to adjust iupseg, nseg, and outseg
                record.nseg = self.renumbering[record.nseg]
                record.iupseg = self.renumbering[record.iupseg]
                record.outseg = self.renumbering[record.outseg]
                ra[ix] = record

            spd[kper] = ra

        self.sfr.segment_data = spd

        # renumber geometry and flow data
        sfr_geometry = self.sfr.channel_geometry_data
        sfr_flow = self.sfr.channel_flow_data

        geometry = {}
        for per, d in sfr_geometry.items():
            t = {}
            for seg, data in d.items():
                t[self.renumbering[seg]] = data

            if len(t) > 0:
                geometry[per] = t

        self.sfr.channel_geometry_data = geometry

        flow = {}
        for per, d in sfr_flow.items():
            t = {}
            for seg, data in d.items():
                t[self.renumbering[seg]] = data

            if len(t) > 0:
                flow[per] = t

        self.sfr.channel_flow_data = flow

    def __renumber_gage(self):
        """
        Reusable method to renumber the GAGE
        package
        """
        if self.gage is None:
            return

        numgage = self.gage.numgage
        gage_data = self.gage.gage_data
        ra = np.recarray((numgage,), dtype=gage_data.dtype)

        for ix, record in enumerate(gage_data):
            if record.gageloc < 0:
                ra[ix] = record
            else:
                record.gageloc = self.renumbering[record.gageloc]
                ra[ix] = record

        self.gage.gage_data = ra

    def __renumber_ag(self):
        """
        Method to renumber the GSFLOW/modflow-nwt
        Agricultural options package
        """
        if self.ag is None:
            return

        if self.ag.irrdiversion is not None:
            irr_data = self.ag.irrdiversion
            spd = {}

            for kper, mflist in irr_data.items():
                ndiv = len(mflist)
                ra = np.recarray((ndiv,), dtype=mflist.dtype)

                for ix, record in enumerate(mflist):
                    record.segid = self.renumbering[record.segid]
                    ra[ix] = record

                spd[kper] = ra

            self.ag.irrdiversion = spd

        if self.ag.supwell is not None:
            sup_data = self.ag.supwell
            spd = {}

            for kper, mflist in sup_data.items():
                nsup = len(mflist)
                ra = np.recarray((nsup,), dtype=mflist.dtype)

                for ix, record in enumerate(mflist):
                    numcell = record.numcell
                    for cell in range(numcell):
                        name = "segid{}".format(cell)
                        record[name] = self.renumbering[record[name]]

                    ra[ix] = record

                spd[kper] = ra

            self.ag.supwell = spd


class Topology(object):
    """
    A topological sort method that uses a modified
    Khan algorithm to sort the SFR network

    Parameters
    ----------
    n_segments : int
        number of sfr segments in network

    """

    def __init__(self, nss=None):
        self.topology = dict()
        self.nss = nss

    def add_connection(self, iseg, ioutseg):
        """
        Method to add a topological connection

        Parmeters
        ---------
        iseg : int
            current segment number
        ioutseg : int
            output segment number
        """
        self.topology[iseg] = ioutseg

    def _sort_util(self, seg, visited, stack):
        """
        Recursive function used by topological
        sort to perform sorting

        Parameters
        ----------
        seg : int
            segment number
        visited : list
            list of bools to indicate if location visited
        stack : list
            stack of sorted segment numbers

        """
        visited[seg] = True
        if seg == 0:
            ioutseg = 0
        else:
            ioutseg = self.topology[seg]

        if not visited[ioutseg]:
            self._sort_util(ioutseg, visited, stack)

        if seg == 0:
            pass
        elif ioutseg == 0:
            stack.append(seg)
        else:
            stack.insert(0, seg)

    def sort(self):
        """
        Method to perform a topological sort
        on the streamflow network

        Returns
        -------
            stack: list of ordered nodes

        """
        visited = {0: False}
        for key in self.topology:
            visited[key] = False

        stack = []
        for i in sorted(visited):
            if i == 0:
                pass
            else:
                if not visited[i]:
                    self._sort_util(i, visited, stack)

        return stack

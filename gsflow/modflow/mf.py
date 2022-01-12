import os
import inspect
import flopy
import platform
from flopy.utils import TemporalReference
from flopy.utils.reference import SpatialReference
from flopy.modflow.mf import Modflow as fpModflow
from ..utils import mfreadnam


class Modflow(fpModflow):
    """
    GSFLOW overridden MODFLOW Model Class. This class is a child of
    the FloPy flopy.modflow.Modflow class

    Parameters
    ----------
    modelname : string, optional
        Name of model.  This string will be used to name the MODFLOW input
        that are created with write_model. (the default is 'modflowtest')
    namefile_ext : string, optional
        Extension for the namefile (the default is 'nam')
    version : string, optional
        Version of MODFLOW to use (the default is 'mf2005').
    exe_name : string, optional
        The name of the executable to use (the default is
        'mf2005').
    listunit : integer, optional
        Unit number for the list file (the default is 2).
    model_ws : string, optional
        model workspace.  Directory name to create model data sets.
        (default is the present working directory).
    external_path : string
        Location for external files (default is None).
    verbose : boolean, optional
        Print additional information to the screen (default is False).
    load : boolean, optional
         (default is True).
    silent : integer
        (default is 0)

    Attributes
    ----------

    Methods
    -------

    See Also
    --------

    Notes
    -----

    Examples
    --------

    >>> import flopy
    >>> m = flopy.modflow.Modflow()

    """

    def __init__(
        self,
        modelname="modflowtest",
        namefile_ext="nam",
        version="mfnwt",
        version2="gsflow",
        exe_name="mfnwt.exe",
        structured=True,
        listunit=2,
        model_ws=".",
        external_path=None,
        verbose=False,
        **kwargs
    ):

        if platform.system().lower() == "windows":
            if not exe_name.endswith(".exe"):
                exe_name += ".exe"

        else:
            if exe_name.endswith(".exe"):
                exe_name = exe_name[:-4]

        super(Modflow, self).__init__(
            modelname=modelname,
            namefile_ext=namefile_ext,
            version=version,
            exe_name=exe_name,
            structured=structured,
            listunit=listunit,
            model_ws=model_ws,
            external_path=external_path,
            verbose=verbose,
            **kwargs
        )

        self.version2 = version2

        # Using a local gsflow import to avoid circular imports with python 2.7
        import gsflow

        # we can override packages here by setting the package object
        # a custom written package if flopy doesn't work with pygsflow
        self.mfnam_packages = {
            "zone": flopy.modflow.ModflowZon,
            "mult": flopy.modflow.ModflowMlt,
            "pval": flopy.modflow.ModflowPval,
            "bas6": flopy.modflow.ModflowBas,
            "dis": flopy.modflow.ModflowDis,
            "disu": flopy.modflow.ModflowDisU,
            "bcf6": flopy.modflow.ModflowBcf,
            "lpf": flopy.modflow.ModflowLpf,
            "hfb6": flopy.modflow.ModflowHfb,
            "chd": flopy.modflow.ModflowChd,
            "fhb": flopy.modflow.ModflowFhb,
            "wel": flopy.modflow.ModflowWel,
            "mnw1": flopy.modflow.ModflowMnw1,
            "mnw2": flopy.modflow.ModflowMnw2,
            "mnwi": flopy.modflow.ModflowMnwi,
            "drn": flopy.modflow.ModflowDrn,
            "rch": flopy.modflow.ModflowRch,
            "evt": flopy.modflow.ModflowEvt,
            "ghb": flopy.modflow.ModflowGhb,
            "gmg": flopy.modflow.ModflowGmg,
            "lmt6": flopy.modflow.ModflowLmt,
            "lmt7": flopy.modflow.ModflowLmt,
            "riv": flopy.modflow.ModflowRiv,
            "str": flopy.modflow.ModflowStr,
            "swi2": flopy.modflow.ModflowSwi2,
            "pcg": flopy.modflow.ModflowPcg,
            "pcgn": flopy.modflow.ModflowPcgn,
            "nwt": flopy.modflow.ModflowNwt,
            "pks": flopy.modflow.ModflowPks,
            "sms": flopy.modflow.ModflowSms,
            "sfr": flopy.modflow.ModflowSfr2,
            "lak": flopy.modflow.ModflowLak,
            "gage": flopy.modflow.ModflowGage,
            "sip": flopy.modflow.ModflowSip,
            "sor": flopy.modflow.ModflowSor,
            "de4": flopy.modflow.ModflowDe4,
            "oc": flopy.modflow.ModflowOc,
            "uzf": flopy.modflow.ModflowUzf1,
            "upw": flopy.modflow.ModflowUpw,
            "sub": flopy.modflow.ModflowSub,
            "swt": flopy.modflow.ModflowSwt,
            "hyd": flopy.modflow.ModflowHyd,
            "hob": flopy.modflow.ModflowHob,
            "vdf": flopy.seawat.SeawatVdf,
            "vsc": flopy.seawat.SeawatVsc,
            "ag": gsflow.modflow.ModflowAg,
        }

    @property
    def model_ws(self):
        return super(Modflow, self).model_ws

    def change_model_ws(self, new_pth=None, reset_external=False):
        """
        Change the model work space.

        Parameters
        ----------
        new_pth : str
            Location of new model workspace.  If this path does not exist,
            it will be created. (default is None, which will be assigned to
            the present working directory).

        Returns
        -------
        val : list of strings
            Can be used to see what packages are in the model, and can then
            be used with get_package to pull out individual packages.

        """
        super(Modflow, self).change_model_ws(new_pth, reset_external)

    def _set_relative_paths(self):
        """
        Fix Gsflow relative paths and names

        """

        if os.path.dirname(self.lst.fn_path) != self.model_ws:
            tmp = self._check_basenames(self.lst)
            self.lst.fn_path = os.path.join(self.model_ws, tmp)
            self.lst.file_name[0] = tmp

        else:
            tmp = self._check_basenames(self.lst)
            self.lst.file_name[0] = tmp
        # rpth = os.path.relpath(self.lst.file_name[0], self.model_ws)
        # self.lst.file_name[0] = rpth

        for pkg in self.packagelist:
            if os.path.dirname(pkg.fn_path) != self.model_ws:
                tmp = self._check_basenames(pkg)
                pkg.fn_path = os.path.join(self.model_ws, tmp)
                pkg.file_name[0] = tmp

            else:
                tmp = self._check_basenames(pkg)
                pkg.file_name[0] = tmp
                pkg.fn_path = os.path.join(self.model_ws, tmp)
            # rpth = os.path.relpath(pkg.file_name[0], self.model_ws)
            # pkg.file_name[0] = rpth

        # skip external fnames for now and then update later
        # for name in self.external_fnames:
        #     rpth = os.path.relpath()

    def _check_basenames(self, pkg):
        """
        Fix to basename not transfering into the fn_path issue

        Parameters
        ----------
        pkg : flopy.modflow.Package

        Returns
        -------
            tmp : str
                new basename for the package
        """
        name = os.path.split(self.name)[-1]
        tmp = os.path.split(pkg.fn_path)[-1]
        if tmp != name + "." + pkg.extension[0]:
            tmp = "{}.{}".format(name, pkg.extension[0])
        return tmp

    def write_name_file(self):
        """
        Write the model name file

        """
        self._set_relative_paths()
        super(Modflow, self).write_name_file()

    def write_input(self, SellPackList=False, check=False):
        """
        Write the input.

        Parameters
        ----------
        SelPackList : False or list of packages

        """
        self._set_relative_paths()
        super(Modflow, self).write_input(SelPackList=SellPackList, check=check)

    @staticmethod
    def load(
        f,
        version="mfnwt",
        exe_name="mfnwt.exe",
        verbose=False,
        model_ws=".",
        load_only=None,
        forgive=False,
        check=True,
        control_file=None,
    ):
        """
        Load an existing MODFLOW model.

        Parameters
        ----------
        f : str
            Path to MODFLOW name file to load.
        version : str, optional
            MODFLOW version. Default 'mf2005', although can be modified on
            loading packages unique to different MODFLOW versions.
        exe_name : str, optional
            MODFLOW executable name. Default 'mf2005.exe'.
        verbose : bool, optional
            Show messages that can be useful for debugging. Default False.
        model_ws : str
            Model workspace path. Default '.' or current directory.
        load_only : list, str or None
            List of case insensitive filetypes to load, e.g. ["bas6", "lpf"].
            One package can also be specified, e.g. "rch". Default is None,
            which attempts to load all files. An empty list [] will not load
            any additional packages than is necessary. At a minimum, "dis" or
            "disu" is always loaded.
        forgive : bool, optional
            Option to raise exceptions on package load failure, which can be
            useful for debugging. Default False.
        check : boolean, optional
            Check model input for common errors. Default True.
        control_file : str, optional
            For using GSFLOW, providing a control file helps to adjust the paths on
            the fly

        Returns
        -------
        ml : Modflow object

        Examples
        --------

        >>> import flopy
        >>> ml = flopy.modflow.Modflow.load('model.nam')

        """

        # similar to modflow command: if file does not exist , try file.nam
        namefile_path = os.path.join(model_ws, f)
        if not os.path.isfile(namefile_path) and os.path.isfile(
            namefile_path + ".nam"
        ):
            namefile_path += ".nam"
        if not os.path.isfile(namefile_path):
            raise IOError("cannot find name file: " + str(namefile_path))

        # Determine model name from 'f', without any extension or path
        modelname = os.path.splitext(os.path.basename(f))[0]

        # if model_ws is None:
        #    model_ws = os.path.dirname(f)
        if verbose:
            print(
                "\nCreating new model with name: {}\n{}\n".format(
                    modelname, 50 * "-"
                )
            )
        ml = Modflow(
            modelname,
            version=version,
            exe_name=exe_name,
            verbose=verbose,
            model_ws=model_ws,
        )

        files_successfully_loaded = []
        files_not_loaded = []

        # set the reference information
        attributes = mfreadnam.attribs_from_namfile_header(namefile_path)

        # read name file
        ext_unit_dict = mfreadnam.parsenamefile(
            namefile_path,
            ml.mfnam_packages,
            control_file=control_file,
            verbose=verbose,
            model_ws=model_ws,
        )
        if ml.verbose:
            print(
                "\n{}\nExternal unit dictionary:\n{}\n{}\n".format(
                    50 * "-", ext_unit_dict, 50 * "-"
                )
            )

        # create a dict where key is the package name, value is unitnumber
        ext_pkg_d = {v.filetype: k for (k, v) in ext_unit_dict.items()}

        # reset version based on packages in the name file
        if "NWT" in ext_pkg_d or "UPW" in ext_pkg_d:
            version = "mfnwt"
        if "GLOBAL" in ext_pkg_d:
            version = "mf2k"
        if "SMS" in ext_pkg_d:
            version = "mfusg"
        if "DISU" in ext_pkg_d:
            version = "mfusg"
            ml.structured = False
        # update the modflow version
        ml.set_version(version)

        # reset unit number for glo file
        if version == "mf2k":
            if "GLOBAL" in ext_pkg_d:
                unitnumber = ext_pkg_d["GLOBAL"]
                filepth = os.path.basename(ext_unit_dict[unitnumber].filename)
                ml.glo.unit_number = [unitnumber]
                ml.glo.file_name = [filepth]
            else:
                ml.glo.unit_number = [0]
                ml.glo.file_name = [""]

        # reset unit number for list file
        if "LIST" in ext_pkg_d:
            unitnumber = ext_pkg_d["LIST"]
            filepth = os.path.basename(ext_unit_dict[unitnumber].filename)
            ml.lst.unit_number = [unitnumber]
            ml.lst.file_name = [filepth]

        # look for the free format flag in bas6
        bas_key = ext_pkg_d.get("BAS6")
        if bas_key is not None:
            bas = ext_unit_dict[bas_key]
            start = bas.filehandle.tell()
            line = bas.filehandle.readline()
            while line.startswith("#"):
                line = bas.filehandle.readline()
            if "FREE" in line.upper():
                ml.free_format_input = True
            bas.filehandle.seek(start)
        if verbose:
            print("ModflowBas6 free format:{0}\n".format(ml.free_format_input))

        # load dis
        dis_key = ext_pkg_d.get("DIS") or ext_pkg_d.get("DISU")
        if dis_key is None:
            raise KeyError("discretization entry not found in nam file")
        disnamdata = ext_unit_dict[dis_key]
        dis = disnamdata.package.load(
            disnamdata.filename, ml, ext_unit_dict=ext_unit_dict, check=False
        )
        files_successfully_loaded.append(disnamdata.filename)
        if ml.verbose:
            print("   {:4s} package load...success".format(dis.name[0]))
        assert ml.pop_key_list.pop() == dis_key
        ext_unit_dict.pop(dis_key)
        start_datetime = attributes.pop("start_datetime", "01-01-1970")
        itmuni = attributes.pop("itmuni", 4)
        ref_source = attributes.pop("source", "defaults")

        if ref_source != "usgs.model.reference":
            itmuni = dis.itmuni
            attributes["lenuni"] = dis.lenuni

        dis.tr = TemporalReference(
            itmuni=itmuni, start_datetime=start_datetime
        )
        dis.start_datetime = start_datetime

        if load_only is None:
            # load all packages/files
            load_only = ext_pkg_d.keys()
        else:  # check items in list
            if not isinstance(load_only, list):
                load_only = [load_only]
            not_found = []
            for i, filetype in enumerate(load_only):
                load_only[i] = filetype = filetype.upper()
                if filetype not in ext_pkg_d:
                    not_found.append(filetype)
            if not_found:
                raise KeyError(
                    "the following load_only entries were not found "
                    "in the ext_unit_dict: " + str(not_found)
                )

        # zone, mult, pval
        if "PVAL" in ext_pkg_d:
            ml.mfpar.set_pval(ml, ext_unit_dict)
            assert ml.pop_key_list.pop() == ext_pkg_d.get("PVAL")
        if "ZONE" in ext_pkg_d:
            ml.mfpar.set_zone(ml, ext_unit_dict)
            assert ml.pop_key_list.pop() == ext_pkg_d.get("ZONE")
        if "MULT" in ext_pkg_d:
            ml.mfpar.set_mult(ml, ext_unit_dict)
            assert ml.pop_key_list.pop() == ext_pkg_d.get("MULT")

        # try loading packages in ext_unit_dict
        for key, item in ext_unit_dict.items():
            if item.package is not None:
                if item.filetype in load_only:
                    if forgive:
                        try:
                            package_load_args = list(
                                inspect.getfullargspec(item.package.load)
                            )[0]
                            if "check" in package_load_args:
                                pck = item.package.load(
                                    item.filename,
                                    ml,
                                    ext_unit_dict=ext_unit_dict,
                                    check=False,
                                )
                            else:
                                pck = item.package.load(
                                    item.filename,
                                    ml,
                                    ext_unit_dict=ext_unit_dict,
                                )
                            files_successfully_loaded.append(item.filename)
                            if ml.verbose:
                                print(
                                    "   {:4s} package load...success".format(
                                        item.filetype
                                    )
                                )
                        except Exception as e:
                            ml.load_fail = True
                            if ml.verbose:
                                print(
                                    "{:4s} package load...failed\n   {!s}".format(
                                        item.filetype, e
                                    )
                                )
                            files_not_loaded.append(item.filename)
                    else:
                        package_load_args = list(
                            inspect.getfullargspec(item.package.load)
                        )[0]
                        if "check" in package_load_args:
                            pck = item.package.load(
                                item.filename,
                                ml,
                                ext_unit_dict=ext_unit_dict,
                                check=False,
                            )
                        else:
                            pck = item.package.load(
                                item.filename, ml, ext_unit_dict=ext_unit_dict
                            )
                        files_successfully_loaded.append(item.filename)
                        if ml.verbose:
                            print(
                                "   {:4s} package load...success".format(
                                    item.filetype
                                )
                            )
                else:
                    if ml.verbose:
                        print(
                            "   {:4s} package load...skipped".format(
                                item.filetype
                            )
                        )
                    files_not_loaded.append(item.filename)
            elif "data" not in item.filetype.lower():
                files_not_loaded.append(item.filename)
                if ml.verbose:
                    print(
                        "   {:4s} package load...skipped".format(item.filetype)
                    )
            elif "data" in item.filetype.lower():
                if ml.verbose:
                    print(
                        "   {} file load...skipped\n      {}".format(
                            item.filetype, os.path.basename(item.filename)
                        )
                    )
                if key not in ml.pop_key_list:
                    # do not add unit number (key) if it already exists
                    if key not in ml.external_units:
                        ml.external_fnames.append(item.filename)
                        ml.external_units.append(key)
                        ml.external_binflag.append(
                            "binary" in item.filetype.lower()
                        )
                        ml.external_output.append(False)
            else:
                raise KeyError("unhandled case: {}, {}".format(key, item))

        # pop binary output keys and any external file units that are now
        # internal
        for key in ml.pop_key_list:
            try:
                ml.remove_external(unit=key)
                ext_unit_dict.pop(key)
            except KeyError:
                if ml.verbose:
                    print(
                        "Warning: external file unit {} does not exist in "
                        "ext_unit_dict.".format(key)
                    )

        # write message indicating packages that were successfully loaded
        if ml.verbose:
            print("")
            print(
                "The following {0} packages were successfully loaded.".format(
                    len(files_successfully_loaded)
                )
            )
            for fname in files_successfully_loaded:
                print("      " + os.path.basename(fname))
            if len(files_not_loaded) > 0:
                print(
                    "   The following {0} packages were not loaded.".format(
                        len(files_not_loaded)
                    )
                )
                for fname in files_not_loaded:
                    print("      " + os.path.basename(fname))
        if check:
            ml.check(f="{}.chk".format(ml.name), verbose=ml.verbose, level=0)

        # return model object
        return ml

# -*- coding: utf-8 -*-
import os
from .control import ControlFile
from .prms import PrmsModel
from .utils import gsflow_io, GsConstant
from .prms import Helper
from .modflow import Modflow
from .modsim import Modsim
import flopy
import subprocess as sp
import platform
import warnings

warnings.simplefilter("always", PendingDeprecationWarning)
warnings.simplefilter("always", UserWarning)


class GsflowModel(object):
    """
    GsflowModel is the GSFLOW model object. This class can be used
    to build a GSFLOW model, to load a GSFLOW model from it's control file,
    to write input files for GSFLOW and to run GSFLOW.

    Parameters
    ----------
    control_file : str
        control file path and name
    prms : PrmsModel object
        gsflow.prms.PrmsModel
    mf : Modflow object
        gsflow.modflow.Modflow
    modflow_only : bool
        flag that indicates only Modflow model
    prms_only : bool
        flag that indicates only PRMS model
    gsflow_exe : str
        GSFLOW executable path and name
    modsim : bool
        boolean flag to indicate that modsim is active
        this creates a gsflow.modsim.Modsim object
    model_ws : str, None
        override method to set the base model directory when the
        GSFLOW control file is not located in the same directory as
        the script to run GSFLOW

    Examples
    --------

    load from control file

    >>> import gsflow
    >>> gsf = gsflow.GsflowModel.load_from_file("gsflow.control")

    create new, empty gsflow object

    >>> control = gsflow.ControlFile(records_list=[])
    >>> gsf = gsflow.GsflowModel(control=control)

    """

    def __init__(
        self,
        control=None,
        prms=None,
        mf=None,
        modflow_only=False,
        prms_only=False,
        gsflow_exe=None,
        modsim=False,
        model_ws=None,
    ):

        if not isinstance(control, ControlFile):
            raise ValueError("control must be a ControlFile object")

        self.control = control
        self.control_file = os.path.abspath(control.control_file)
        self.ws = None
        self._modflow_only = modflow_only
        self._prms_only = prms_only
        self.prms = None
        self.mf = None
        self.modsim = None
        self.gsflow_exe = gsflow_exe

        if gsflow_exe is None:
            self.gsflow_exe = os.path.join(
                os.path.dirname(__file__), r"bin\gsflow.exe"
            )

        # set prms object
        if not modflow_only:
            if prms and isinstance(prms, PrmsModel):
                self.prms = prms
            else:
                err = "prms is not a PrmsModel object, skipping..."
                warnings.warn(err, UserWarning)

        # set flopy modflow object
        if not prms_only:
            if mf and isinstance(mf, flopy.modflow.Modflow):
                self.mf = mf
                namefile = os.path.basename(
                    control.get_values("modflow_name")[0]
                )
                if namefile is not None:
                    self.mf.namefile = namefile
            else:
                err = "modflow is not a gsflow.modflow.Modflow object, skipping..."
                warnings.warn(err, UserWarning)

        if modsim:
            self.modsim = Modsim(self)

        self.help = Helper()

    @property
    def modflow_only(self):
        """
        Returns
        -------
            bool
        """
        return self._modflow_only

    @property
    def prms_only(self):
        """
        Returns
        -------
            bool
        """
        return self._prms_only

    def export_nc(self, f, **kwargs):
        """
        Method to export the GSFLOW model as a netcdf
        file. This method only works if nhru is equivalent
        to nrow * ncol in modflow.

        Parameters
        ----------
        f : str
            netcdf file name
        kwargs :
            keyword arguments for netcdf

        """
        if not f.endswith(".nc"):
            raise AssertionError("f must end with .nc extension")
        if self.mf is None:
            err = "Modflow object must be loaded to export netcdf file"
            raise AssertionError(err)

        f = self.mf.export(f, **kwargs)
        if self.prms is not None:
            f = self.prms.export_nc(f, self.mf, **kwargs)

        return f

    @staticmethod
    def load_from_file(
        control_file,
        gsflow_exe="gsflow.exe",
        modflow_only=False,
        prms_only=False,
        mf_load_only=None,
        forgive=False,
        model_ws=None,
    ):
        """
        Method to load a gsflow model from it's control file

        Parameters
        ----------
        control_file : str
            control file path & name, GSFLOW
        gsflow_exe : str
            gsflow executable path & name
        modflow_only : bool
            flag to load only modflow from the control file
        prms_only : bool
            flag to load only prms from the control file
        mf_load_only : list
            list of packages to load from modflow ex. [DIS, BAS, LPF]
        forgive : bool
            forgive file loading errors in flopy
        model_ws : str, None
            override method to set the base model directory when the
            GSFLOW control file is not located in the same directory as
            the script to run GSFLOW
        Returns
        -------
            GsflowModel object

        Examples
        --------

        >>> import gsflow
        >>> gsf = gsflow.GsflowModel.load_from_file("gsflow.control")

        """
        prms = None
        modflow = None
        modsim = False
        if not (os.path.isfile(control_file)):
            raise ValueError("Cannot find control file")

        if model_ws is not None:
            control = ControlFile.load_from_file(control_file, abs_path=False)
        else:
            control = ControlFile.load_from_file(control_file)
        print("Control file is loaded")

        mode = control.get_values("model_mode")[0].upper()
        if mode == "MODFLOW":
            modflow_only = True
        elif mode == "PRMS":
            prms_only = True
        elif "MODSIM" in mode:
            modsim = True
        else:
            pass

        # load prms
        if not modflow_only:
            print("Working on loading PRMS model ...")
            prms = PrmsModel.load_from_file(control_file, model_ws=model_ws)

        if not prms_only:
            # get model mode
            if "GSFLOW" in mode.upper() or "MODFLOW" in mode.upper():
                print("Working on loading MODFLOW files ....")
                modflow = GsflowModel._load_modflow(
                    control, mf_load_only, model_ws, forgive
                )
                print("MODFLOW files are loaded ... ")

            else:
                prms_only = True
                modflow_only = False
                print("Mode is set to PRMS only, loading PRMS model only")

        return GsflowModel(
            control=control,
            prms=prms,
            mf=modflow,
            modflow_only=modflow_only,
            prms_only=prms_only,
            gsflow_exe=gsflow_exe,
            modsim=modsim,
        )

    @staticmethod
    def _load_modflow(control, mf_load_only, model_ws=None, forgive=False):
        """
        The package files in the .nam file are relative to the execuatble
        gsflow. Here we set the model_ws to the location of the gsflow exe, via
        the control file or a user supplied model_ws parameter

        Parameters
        ----------
        control : ControlFile object
            control file object
        mf_load_only : list
            list of packages to restrict modflow loading to
        model_ws : str
            optional parameter that allows the use to set the model_ws
        forgive : bool
            forgive file load errors in modflow


        Returns
        -------
            Modflow object

        """
        name = control.get_values("modflow_name")
        control_file = control.control_file
        if model_ws is None:
            name = gsflow_io.get_file_abs(
                control_file=control_file, fn=name[0]
            )
            model_ws, name = os.path.split(name)
        else:
            model_ws = gsflow_io.get_file_abs(model_ws=model_ws)
            name = name[0]
            control_file = None

        return Modflow.load(
            name,
            model_ws=model_ws,
            control_file=control_file,
            load_only=mf_load_only,
            forgive=forgive,
        )

    def write_input(self, basename=None, workspace=None, write_only=None):
        """
         Write input files for gsflow. Four cases are possible:
            (1) if basename and workspace are None,then the exisiting files will be overwritten
            (2) if basename is specified, only file names will be changes
            (3) if only workspace is specified, only folder will be changed
            (4) when both basename and workspace are specifed both files are changed

        Parameters
        ----------
        basename : str
            project basename
        workspace :  str
            model output directory
        write_only: a list
            ['control', 'parameters', 'prms_data', 'mf', 'modsim']

        Examples
        --------

        >>> gsf = gsflow.GsflowModel.load_from_file('gsflow.control')
        >>> gsf.write_input(basename="new", workspace="../new_model")

        """
        print("Writing the project files .....")
        if workspace is not None:
            workspace = os.path.abspath(workspace)

        if (basename, workspace) == (None, None):
            print("Warning: input files will be overwritten....")
            self._write_all(write_only)

        # only change the directory
        elif basename is None and workspace is not None:
            if not (os.path.isdir(workspace)):
                os.mkdir(workspace)
            fnn = os.path.basename(self.control.control_file)
            self.control.model_dir = workspace
            self.control.control_file = os.path.join(workspace, fnn)
            self.control_file = os.path.join(workspace, fnn)
            if self.prms is not None:
                self.prms.control_file = self.control_file

                # change parameters
                new_param_file_list = []
                for par_record in self.prms.parameters.parameters_list:
                    curr_file = os.path.basename(par_record.file_name)
                    curr_file = os.path.join(workspace, curr_file)
                    par_record.file_name = curr_file
                    if not (curr_file in new_param_file_list):
                        new_param_file_list.append(curr_file)
                self.control.set_values("param_file", new_param_file_list)

                # change datafile
                curr_file = os.path.relpath(
                    os.path.join(workspace, self.prms.data.name),
                    self.control.model_dir,
                )
                self.prms.data.model_dir = workspace
                self.control.set_values("data_file", [curr_file])

            # change mf
            if self.mf is not None:
                self.mf.change_model_ws(workspace, reset_external=True)
                mfnm = self.mf.name + ".nam"
                self.control.set_values("modflow_name", [mfnm])

            # update file names in control object
            self._update_control_fnames(workspace, basename)
            # write
            if self.prms is not None:
                self.prms.control = self.control
            self._write_all(write_only)

        # only change the basename
        elif basename is not None and workspace is None:
            cnt_file = basename + "_cont.control"
            ws_ = os.path.dirname(self.control.control_file)
            self.control.control_file = os.path.join(ws_, cnt_file)
            self.control_file = os.path.join(ws_, cnt_file)
            self.prms.control_file = self.control_file

            # change parameters
            flist = self.prms.parameters.parameter_files
            new_param_file_list = []
            for ifile, par_record in enumerate(
                self.prms.parameters.parameters_list
            ):
                file_index = flist.index(par_record.file_name)
                par_file = basename + "_par_{}.params".format(file_index)
                curr_dir = self.control.model_dir
                curr_file = os.path.join(curr_dir, par_file)
                par_record.file_name = curr_file
                if not (curr_file in new_param_file_list):
                    new_param_file_list.append(curr_file)
            self.control.set_values("param_file", new_param_file_list)

            # change datafile
            dfile = basename + "_dat.data"
            curr_file = os.path.relpath(
                os.path.join(self.prms.data.model_dir, dfile),
                self.control.model_dir,
            )
            self.prms.data.name = dfile
            self.control.set_values("data_file", [curr_file])

            # change mf
            if self.mf is not None:
                curr_dir = self.mf.model_ws
                self.mf._set_name(basename)
                self._update_mf_basename(basename)
                mfnm = self.mf.name + ".nam"
                self.control.set_values("modflow_name", [mfnm])

            # update file names in control object
            self._update_control_fnames(workspace, basename)
            self.prms.control = self.control
            self._write_all(write_only)

        # change both directory & basename
        elif basename is not None and workspace is not None:
            if not (os.path.isdir(workspace)):
                os.mkdir(workspace)
            cnt_file = basename + "_cont.control"
            self.control.model_dir = workspace
            self.control.control_file = os.path.join(workspace, cnt_file)
            self.prms.control_file = self.control.control_file
            self.control_file = self.control.control_file

            # change parameters
            # get param files list
            flist = self.prms.parameters.parameter_files
            new_param_file_list = []
            for ifile, par_record in enumerate(
                self.prms.parameters.parameters_list
            ):
                file_index = flist.index(par_record.file_name)
                par_file = basename + "_par_{}.params".format(file_index)
                curr_file = os.path.join(workspace, par_file)
                par_record.file_name = curr_file
                if not (curr_file in new_param_file_list):
                    new_param_file_list.append(curr_file)
            self.control.set_values("param_file", new_param_file_list)
            # change datafile
            dfile = basename + "_dat.data"
            curr_file = os.path.relpath(
                os.path.join(workspace, dfile), self.control.model_dir
            )
            self.prms.data.model_dir = workspace
            self.prms.data.name = dfile
            self.control.set_values("data_file", [curr_file])

            # flatten mf
            if self.mf is not None:
                self.mf.change_model_ws(workspace)
                self.mf._set_name(os.path.join(workspace, basename))
                self._update_mf_basename(basename)

            mfnm = basename + ".nam"
            self.control.set_values(
                "modflow_name",
                [
                    os.path.relpath(
                        os.path.join(workspace, mfnm), self.control.model_dir
                    )
                ],
            )
            # update file names in control object
            self._update_control_fnames(workspace, basename)
            self.prms.control = self.control
            self._write_all(write_only)

        else:
            raise NotImplementedError()

    def _update_control_fnames(self, workspace, basename):
        """
        Method to update control file names and paths

        Parameters
        ----------
        workspace :  str
            model output directory
        basename : str
            project basename

        """
        if workspace is not None and basename is None:
            self.control.model_dir = workspace
            for rec_name in GsConstant.GSFLOW_FILES:
                if rec_name in self.control.record_names:
                    file_values = self.control.get_values(rec_name)
                    file_value = []
                    for fil in file_values:
                        va = os.path.join(workspace, os.path.basename(fil))
                        va = os.path.relpath(va, self.control.model_dir)
                        file_value.append(va)
                    self.control.set_values(rec_name, file_value)

        else:
            for rec_name in GsConstant.GSFLOW_FILES:
                if rec_name in self.control.record_names:
                    if rec_name in ("modflow_name",):
                        continue

                    elif rec_name in (
                        "modflow_name",
                        "param_file",
                        "data_file",
                    ):
                        file_values = self.control.get_values(rec_name)
                        file_value = []
                        for fil in file_values:
                            ws, filvalue = os.path.split(fil)
                            if not ws:
                                pass
                            else:
                                filvalue = os.path.relpath(
                                    fil, self.control.model_dir
                                )

                            file_value.append(filvalue)
                        self.control.set_values(rec_name, file_value)

                    else:
                        file_values = self.control.get_values(rec_name)
                        file_value = []
                        for fil in file_values:
                            if workspace is None:
                                workspace = self.control.model_dir
                            vvfile = rec_name.split("_")
                            del vvfile[-1]
                            vvfile = "_".join(vvfile)
                            if "." in fil:
                                ext = fil.split(".")[-1]
                            else:
                                ext = "dat"
                            vvfile = basename + "_" + vvfile + "." + ext
                            filvalue = os.path.join(workspace, vvfile)
                            filvalue = os.path.relpath(
                                filvalue, self.control.model_dir
                            )
                            file_value.append(filvalue)
                        self.control.set_values(rec_name, file_value)

    def _update_mf_basename(self, basename):
        """
        Convience method to update modflow Basename

        Parameters
        ----------
        basename : str
            basename of the Modflow object

        """
        out_files_list = []
        for ix, out_file in enumerate(self.mf.output_fnames):

            if out_file.count(".") > 1:
                ext = out_file.split(".")
                del ext[0]
                ext = ".".join(ext)
            else:
                ext = out_file.split(".")[-1]

            new_outfn = "{}.{}".format(basename, ext)
            out_files_list.append(new_outfn)
        self.mf.output_fnames = out_files_list

    def _write_all(self, write_only):
        """
        Method to write input files

        Parameters
        ----------
        write_only : list
            list of files to write accepts,
            control, parameters, prms_data, mf, and modsim

        """

        write_only_options = (
            "control",
            "parameters",
            "prms_data",
            "mf",
            "modsim",
        )
        if write_only is not None:
            if not isinstance(write_only, list):
                raise ValueError("write_only agrgument must be a list")

            # make write options case insensitive
            write_only = [i.lower() for i in write_only]
            for write_option in write_only:
                if not (write_option in write_only_options):
                    raise ValueError(
                        "The option '{}' is not recognized...".format(
                            write_option
                        )
                    )
        else:
            write_only = ()

        # write control
        if len(write_only) == 0 or "control" in write_only:
            print("Writing Control file ...")
            self.control.write()

        if self.prms is not None:
            # self write parameters
            if len(write_only) == 0 or "parameters" in write_only:
                print("Writing Parameters files ...")
                self.prms.parameters.write()

            # write data
            if len(write_only) == 0 or "prms_data" in write_only:
                print("Writing Data file ...")
                self.prms.data.write()

        # write mf
        if self.mf is not None:
            if len(write_only) == 0 or "mf" in write_only:
                print("Writing Modflow files...")
                self.mf.write_input()

        if self.modsim is not None:
            if len(write_only) == 0 or "modsim" in write_only:
                print("Writing MODSIM shapefile")
                self.modsim.write_modsim_shapefile()

    def run_model(self, model_ws=".", forgive=False):
        """
        Method to run a gsflow model

        Parameters
        ----------
        model_ws : str
            parameter to specify the model directory
        forgive : bool
            forgives convergence issues

        Returns
        -------
            None or (success, buffer)

        Examples
        --------

        >>> gsf = gsflow.GsflowModel.load_from_file("gsflow.control")
        >>> gsf.run_model()

        """
        fn = self.control_file
        if not os.path.isfile(self.gsflow_exe):
            print(
                "Warning : The executable of the model is not specified. Use .gsflow_exe "
                "to define its path... "
            )
            return None

        normal_msg = [
            "normal termination",
        ]  # , "simulation successful"]
        if forgive:
            normal_msg.append("failed to meet solver convergence criteria")

        return self.__run(
            exe_name=self.gsflow_exe,
            namefile=fn,
            normal_msg=normal_msg,
            model_ws=model_ws,
        )

    def _generate_batch_file(self):
        fn = os.path.dirname(self.control_file)
        fn = os.path.join(fn, "__run_gsflow.bat")
        self.__bat_file = fn
        fidw = open(fn, "w")
        exe = os.path.normpath(os.path.join(os.getcwd(), self.gsflow_exe))
        cmd = exe + " " + self.control_file
        fidw.write(cmd)
        fidw.close()

    def __run(
        self,
        exe_name,
        namefile,
        model_ws=".",
        silent=False,
        report=False,
        normal_msg="normal termination",
        cargs=None,
    ):
        """
        This function will run the model using subprocess.Popen.

        Parameters
        ----------
        exe_name : str
            Executable name (with path, if necessary) to run.
        namefile : str
            Namefile of model to run. The namefile must be the
            filename of the namefile without the path.
        model_ws : str
            Path to the location of the namefile. (default is the
            current working directory - './')
        silent : boolean
            Echo run information to screen (default is True).
        report : boolean, optional
            Save stdout lines to a list (buff) which is returned
            by the method . (default is False).
        normal_msg : str
            Normal termination message used to determine if the
            run terminated normally. (default is 'normal termination')
        cargs : str or list of strings
            additional command line arguments to pass to the executable.
            Default is None
        Returns
        -------
        (success, buff)
        success : boolean
        buff : list of lines of stdout

        """

        def is_exe(fpath):
            return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

        def which(program):
            fpath, fname = os.path.split(program)
            if fpath:
                if is_exe(program):
                    return program
            else:
                # test for exe in current working directory
                if is_exe(program):
                    return program
                # test for exe in path statement
                for path in os.environ["PATH"].split(os.pathsep):
                    path = path.strip('"')
                    exe_file = os.path.join(path, program)
                    if is_exe(exe_file):
                        return exe_file
            return None

        success = False
        buff = []

        # convert normal_msg to lower case for comparison
        if isinstance(normal_msg, str):
            normal_msg = [normal_msg.lower()]
        elif isinstance(normal_msg, list):
            for idx, s in enumerate(normal_msg):
                normal_msg[idx] = s.lower()

        # Check to make sure that program and namefile exist
        exe = which(exe_name)
        if exe is None:
            if platform.system() in "Windows":
                if not exe_name.lower().endswith(".exe"):
                    exe = which(exe_name + ".exe")
        if exe is None:
            s = "The program {} does not exist or is not executable.".format(
                exe_name
            )
            raise Exception(s)
        else:
            if not silent:
                s = "pyGSFLOW is using the following executable to run the model: {}".format(
                    exe
                )
                print(s)

        exe = os.path.normpath(os.path.join(os.getcwd(), exe))

        if not os.path.isfile(os.path.join(model_ws, namefile)):
            s = "The namefile for this model does not exists: {}".format(
                namefile
            )
            raise Exception(s)

        # simple little function for the thread to target
        #  def q_output(output, q):
        #    for line in iter(output.readline, b''):
        #        q.put(line)
        #        time.sleep(1)
        #        output.close()

        # create a list of arguments to pass to Popen

        argv = [exe, namefile]

        # add additional arguments to Popen arguments
        if cargs is not None:
            if isinstance(cargs, str):
                cargs = [cargs]
            for t in cargs:
                argv.append(t)

        # run the model with Popen
        if platform.system().lower() == "windows":
            self._generate_batch_file()
            argv = self.__bat_file
        else:
            pass

        model_ws = os.path.dirname(self.control_file)
        proc = sp.Popen(argv, stdout=sp.PIPE, stderr=sp.STDOUT, cwd=model_ws)

        while True:
            line = proc.stdout.readline()
            c = line.decode("utf-8")
            if c != "":
                for msg in normal_msg:
                    if msg in c.lower():
                        success = True
                        break
                c = c.rstrip("\r\n")
                if not silent:
                    print("{}".format(c))
                if report:
                    buff.append(c)
            else:
                break
        return success, buff

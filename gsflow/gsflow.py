# -*- coding: utf-8 -*-
import os, sys
import logging
from .control import ControlFile
from .prms import PrmsModel
from .utils import io, GsConstant
from .prms import Helper
import flopy
from .modflow import Modflow
import subprocess as sp
import platform

if sys.version_info > (3, 0):
    import queue as Queue
else:
    import Queue
from datetime import datetime
import threading
import warnings
warnings.simplefilter('always', PendingDeprecationWarning)
warnings.simplefilter('always', UserWarning)


# todo: make this a static method
def load(control_file):
    gs = Gsflow(control_file=control_file)
    gs.load()
    return gs


class Gsflow(object):
    def __new__(cls, control_file=None, prms=None, mf=None, mf_load_only=None,
                 prms_load_only=None, gsflow_exe=None):
        err = "Gsflow has been deprecated, calling GsflowModel.load_from_file()"
        warnings.warn(err, PendingDeprecationWarning)
        return GsflowModel.load_from_file(control_file=control_file, gsflow_exe=gsflow_exe,
                                          modflow_only=mf_load_only, prms_only=prms_load_only)


class GsflowModel(object):
    """

    Parameters
    ----------
    control_file
    prms
    mf
    mf_load_only
    prms_load_only
    gsflow_exe

    """
    def __init__(self, control=None, prms=None, mf=None, modflow_only=False,
                 prms_only=False, gsflow_exe=None):

        if not isinstance(control, ControlFile):
            raise ValueError("control must be a ControlFile object")

        self.control = control
        self.control_file = os.path.abspath(control.control_file)
        self.ws = None
        self._modflow_only = modflow_only
        self._prms_only = prms_only
        self.prms = None
        self.mf = None
        self.gsflow_exe = gsflow_exe

        if gsflow_exe is None:
            self.gsflow_exe = os.path.join(os.path.dirname(__file__), r"bin\gsflow.exe")

        # set prms modflow object
        if not modflow_only:
            if prms and isinstance(prms, PrmsModel):
                self.prms = prms
            else:
                err = "prms is not a PrmsModel object, skipping..."
                warnings.warn(err, UserWarning)

        # set flopy modflow object
        if not prms_only:
            # todo: change this to gsflow.modflow.Modflow instance
            if mf and isinstance(mf, flopy.modflow.Modflow):
                self.mf = mf
                namefile = os.path.basename(control.get_values('modflow_name')[0])
                if namefile is not None:
                    self.mf.namefile = namefile
            else:
                err = "prms is not a gsflow.modflow.Modflow object, skipping..."
                warnings.warn(err, UserWarning)

        self.help = Helper()

    @property
    def modflow_only(self):
        return self._modflow_only

    @property
    def prms_only(self):
        return self._prms_only

    @property
    def mf_load_only(self):
        err = "mf_load_only is deprecated, calling modflow_only"
        warnings.warn(err, PendingDeprecationWarning)
        return self.modflow_only

    @property
    def prms_load_only(self):
        err = "prms_load only is deprecated, calling prms_only"
        warnings.warn(err, PendingDeprecationWarning)
        return self.prms_only

    @property
    def Help(self):
        err = "Help is deprecated, calling help"
        warnings.warn(err, PendingDeprecationWarning)
        return self.help

    @staticmethod
    def load_from_file(control_file, gsflow_exe="gsflow.exe", modflow_only=False,
                       prms_only=False):
        """

        Parameters
        ----------
        control_file
        gsflow_exe
        modflow_only
        prms_only

        Returns
        -------

        """
        prms = None
        modflow = None
        if not (os.path.isfile(control_file)):
            raise ValueError("Cannot find control file")

        control = ControlFile.load_from_file(control_file)
        print("Control file is loaded")

        # load prms
        if not modflow_only:
            print("Working on loading PRMS model ...")
            prms = PrmsModel.load_from_file(control_file)

        if not prms_only:
            # get model mode
            mode = control.get_values('model_mode')
            if 'GSFLOW' in mode[0] or 'MODFLOW' in mode[0]:
                print("Working on loading MODFLOW files ....")
                modflow = GsflowModel._load_modflow(control)
                print("MOSFLOW files are loaded ... ")
                namefile = os.path.basename(control.get_values('modflow_name')[0])
            else:
                prms_only = True
                modflow_only = False
                print("Mode is set to PRMS only, loading PRMS model only")

        return GsflowModel(control=control, prms=prms, mf=modflow,
                           modflow_only=modflow_only, prms_only=prms_only,
                           gsflow_exe=gsflow_exe)

    @staticmethod
    def _load_modflow(control):
        """
        The package files in the .nam file are relative to the execuatble gsflow. So here, we generate a temp.nam
        file that that has the absolute files

        We should be able to deprecate this after some time by overloading the Modflow load scheme and
        the modflow write_input scheme

        Parameters
        ----------
        fname : str
            file name of the modflow name file
        control : ControlFile object
            control file object

        Returns
        -------

        """
        name = control.get_values('modflow_name')
        control_file = control.control_file
        name = io.get_file_abs(control_file=control_file, fn=name[0])
        model_dir, name = os.path.split(name)
        return Modflow.load(name, model_ws=model_dir, control_file=control_file)

    # def change_ws(self, ws):
    #
    #     if os.path.isdir(ws):
    #         print("Warning: The {} directory already exists".format(ws))
    #     parent_folder = os.path.dirname(ws)
    #
    #     if not (os.path.isdir(parent_folder)):
    #         raise ValueError(" The parent directory {} doesn't exist...".format(parent_folder))
    #
    #     if not (os.path.isdir(ws)):
    #         os.mkdir(ws)
    #
    #     self.ws = ws
    #
    #     # change control file location
    #     fnn = os.path.basename(self.control.control_file)
    #     self.control.control_file = os.path.join(self.ws, fnn)
    #
    #     # change parameters
    #     for par_record in self.prms.Parameters.parameters_list:
    #         curr_file = os.path.basename(par_record.file_name)
    #         curr_file = os.path.join(self.ws, curr_file)
    #         par_record.file_name = curr_file
    #
    #     # change datafile
    #     curr_file = os.path.basename(self.prms.Data.data_file)
    #     curr_file = os.path.join(self.ws, curr_file)
    #     self.prms.Data.data_file = curr_file
    #
    #     # change mf
    #     if not (self.mf == None):
    #         self.mf.change_model_ws(self.ws)

    # def change_base_file_name(self, filename):
    #     # change control file location
    #     cnt_file = filename + "_cnt" + ".control"
    #     dir__ = os.path.dirname(self.control.control_file)
    #     self.control.control_file = os.path.join(dir__, cnt_file)
    #
    #     # change parameters
    #     for index, par_record in enumerate(self.prms.Parameters.parameters_list):
    #         curr_file = os.path.basename(par_record.file_name)
    #         curr_file = os.path.join(self.ws, curr_file)
    #         par_record.file_name = curr_file
    #
    #     # change datafile
    #     curr_file = os.path.basename(self.prms.Data.data_file)
    #     curr_file = os.path.join(self.ws, curr_file)
    #     self.prms.Data.data_file = curr_file
    #     pass


    # def _mk_dir(self, dir_):
    #     if not (os.path.isdir(dir_)):
    #         os.mkdir(dir_)
    #     else:
    #         print(" Warning:  the directory exists {}".format(dir_))

    def write_input(self, basename=None, workspace=None):
        """
         Write input files for gsflow. Four cases are possible:
            (1) if basename and workspace are None,then the exisiting files will be overwritten
            (2) if basename is specified, only file names will be changes
            (3) if only workspace is specified, only folder will be changed
            (4) when both basename and workspace are specifed both files are changed

        Parameters
        ----------
        basename
        workspace

        Returns
        -------

        """
        # overwrite

        print("Writing the project files .....")
        if workspace is not None:
            workspace = os.path.abspath(workspace)

        if (basename, workspace) == (None, None):
            print("Warning: input files will be overwritten....")
            self._write_all()
            return

        # only change the directory
        if basename is None and workspace is not None:
            if not (os.path.isdir(workspace)):
                os.mkdir(workspace)
            fnn = os.path.basename(self.control.control_file)
            self.control.model_dir = workspace
            self.control.control_file = os.path.join(workspace, fnn)
            self.control_file = os.path.join(workspace, fnn)
            self.prms.control_file = self.control_file

            # change parameters
            new_param_file_list = []
            for par_record in self.prms.parameters.parameters_list:
                curr_file = os.path.basename(par_record.file_name)
                curr_file = os.path.join(workspace, curr_file)
                par_record.file_name = curr_file
                if not (curr_file in new_param_file_list):
                    new_param_file_list.append(curr_file)
            self.control.set_values('param_file', new_param_file_list)

            # change datafile
            curr_file = os.path.relpath(os.path.join(workspace, self.prms.data.name),
                                        self.control.model_dir)
            self.prms.data.model_dir = workspace
            self.control.set_values('data_file', [curr_file])

            # change mf
            if self.mf is not None:
                self.mf.change_model_ws(workspace)
                mfnm = self.mf.name + ".nam"
                self.control.set_values('modflow_name', [mfnm])

            # update file names in control object
            self._update_control_fnames(workspace, basename)
            # write
            self.prms.control = self.control
            self._write_all()
            return

        # only change the basename
        if basename is not None and workspace is None:
            cnt_file = basename + "_cont.control"
            ws_ = os.path.dirname(self.control.control_file)
            self.control.control_file = os.path.join(ws_, cnt_file)
            self.control_file = os.path.join(ws_, cnt_file)
            self.prms.control_file = self.control_file

            # change parameters
            flist = self.prms.parameters.parameter_files
            new_param_file_list = []
            for ifile, par_record in enumerate(self.prms.parameters.parameters_list):
                file_index = flist.index(par_record.file_name)
                par_file = basename + "_par_{}.params".format(file_index)
                curr_dir = self.control.model_dir # os.path.dirname(par_record.file_name)
                curr_file = os.path.join(curr_dir, par_file)
                par_record.file_name = curr_file
                if not (curr_file in new_param_file_list):
                    new_param_file_list.append(curr_file)
            self.control.set_values('param_file', new_param_file_list)

            # change datafile
            dfile = basename + "_dat.data"
            curr_file = os.path.relpath(os.path.join(self.prms.data.model_dir, dfile),
                                        self.control.model_dir)
            self.prms.data.name = dfile
            self.control.set_values('data_file', [curr_file])

            # change mf
            if self.mf is not None:
                curr_dir = self.mf.model_ws
                self.mf._set_name(basename)
                self._update_mf_basename(basename)
                mfnm = self.mf.name + ".nam"
                self.control.set_values('modflow_name',[mfnm])

            # update file names in control object
            self._update_control_fnames(workspace, basename)
            self.prms.control = self.control
            self._write_all()
            return

        # change both directory & basename
        if basename is not None and workspace is not None:
            if not (os.path.isdir(workspace)):
                os.mkdir(workspace)
            cnt_file = basename + "_cont.control"
            self.control.model_dir = workspace
            self.control.control_file = os.path.join(workspace, cnt_file)
            self.prms.control_file = self.control.control_file
            self.control_file = self.control.control_file

            # change parameters
            ## get param files list
            flist = self.prms.parameters.parameter_files
            new_param_file_list = []
            for ifile, par_record in enumerate(self.prms.parameters.parameters_list):
                file_index = flist.index(par_record.file_name)
                par_file = basename + "_par_{}.params".format(file_index)
                curr_file = os.path.join(workspace, par_file)
                par_record.file_name = curr_file
                if not (curr_file in new_param_file_list):
                    new_param_file_list.append(curr_file)
            self.control.set_values('param_file', new_param_file_list)
            # change datafile
            dfile = basename + "_dat.data"
            curr_file = os.path.relpath(os.path.join(workspace, dfile),
                                        self.control.model_dir)
            self.prms.data.model_dir = workspace
            self.prms.data.name = dfile
            self.control.set_values('data_file', [curr_file])

            # flatten mf
            if self.mf is not None:
                self.mf.change_model_ws(workspace)
                self.mf._set_name(os.path.join(workspace, basename))
                self._update_mf_basename(basename)

            mfnm = basename + ".nam"
            self.control.set_values('modflow_name', [os.path.relpath(os.path.join(workspace, mfnm),
                                                                     self.control.model_dir)])
            # update file names in control object
            self._update_control_fnames(workspace, basename)
            self.prms.control = self.control
            self._write_all()
            return

    def _update_control_fnames(self, workspace, basename):
        """
        Method to update control file names and paths

        Parameters
        ----------
        workspace : model output directory
        basename : project basename

        Returns
        -------

        """
        if workspace is not None and basename is None:
            self.control.model_dir = workspace
            for rec_name in GsConstant.GSFLOW_FILES:
                if rec_name in self.control.record_names:
                    file_values = self.control.get_values(rec_name)
                    file_value = []
                    for fil in file_values:
                        cnt_dir = os.path.dirname(self.control_file)
                        va = os.path.join(workspace, os.path.basename(fil))
                        va = os.path.relpath(va, self.control.model_dir)
                        file_value.append(va)
                    self.control.set_values(rec_name, file_value)

        else:
            for rec_name in GsConstant.GSFLOW_FILES:
                if rec_name in self.control.record_names:
                    if rec_name in ('modflow_name'):
                        continue

                    elif rec_name in ('modflow_name', 'param_file', 'data_file'):
                        file_values = self.control.get_values(rec_name)
                        file_value = []
                        for fil in file_values:
                            ws, filvalue = os.path.split(fil)
                            if not ws:
                                pass
                            else:
                                filvalue = os.path.relpath(fil, self.control.model_dir)

                            file_value.append(filvalue)
                        self.control.set_values(rec_name, file_value)

                    else:
                        file_values = self.control.get_values(rec_name)
                        file_value = []
                        for fil in file_values:
                            if workspace is None:
                                workspace = self.control.model_dir # os.path.dirname(fil)
                            vvfile = rec_name.split("_")
                            del vvfile[-1]
                            vvfile = "_".join(vvfile)
                            if "." in fil:
                                ext = fil.split(".")[-1]
                            else:
                                ext = "dat"
                            vvfile = basename + "_" + vvfile + "." + ext
                            filvalue = os.path.join(workspace, vvfile)
                            filvalue = os.path.relpath(filvalue, self.control.model_dir)
                            file_value.append(filvalue)
                        self.control.set_values(rec_name, file_value)

    def _update_mf_basename(self, basename):
        """
        Convience method to update modflow Basename

        Parameters
        ----------
        basename : str

        """
        out_files_list = []
        for ix, out_file in enumerate(self.mf.output_fnames):

            if out_file.count('.') > 1:
                ext = out_file.split(".")
                del ext[0]
                ext = ".".join(ext)
            else:
                ext = out_file.split(".")[-1]

            new_outfn = "{}.{}".format(basename, ext)
            out_files_list.append(new_outfn)
        self.mf.output_fnames = out_files_list

    def _write_all(self):

        # write control
        self.control.write()
        print("Control file is written...")

        # self write parameters
        self.prms.parameters.write()
        print("Parameters files are written...")

        # write data
        self.prms.data.write()
        print("Data file is written...")

        # write mf
        if self.mf is not None:
            self.mf.write_input()
            print("Modflow files are written...")

    def run_model(self):
        fn = self.control_file
        cnt_folder = os.path.dirname(fn)
        fnm = os.path.abspath(fn)
        if not os.path.isfile(self.gsflow_exe):
            print ("Warning : The executable of the model is not specified. Use .gsflow_exe "
                   "to define its path... ")
            return None
        return self.__run(exe_name=self.gsflow_exe, namefile=fn)

    def _generate_batch_file(self):
        fn = os.path.dirname(self.control_file)
        fn = os.path.join(fn, "__run_gsflow.bat")
        self.__bat_file = fn
        fidw = open(fn, 'w')
        exe = os.path.normpath(os.path.join(os.getcwd(), self.gsflow_exe))
        cmd = exe + " " + self.control_file
        fidw.write(cmd)
        fidw.close()

    def __run(self, exe_name, namefile, model_ws='./',
              silent=False, pause=False, report=False,
              normal_msg='normal termination',
              async_run=False, cargs=None):
        """
        This function will run the model using subprocess.Popen.  It
        communicates with the model's stdout asynchronously and reports
        progress to the screen with timestamps

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
        pause : boolean, optional
            Pause upon completion (default is False).
        report : boolean, optional
            Save stdout lines to a list (buff) which is returned
            by the method . (default is False).
        normal_msg : str
            Normal termination message used to determine if the
            run terminated normally. (default is 'normal termination')
        async_run : boolean
            asynchonously read model stdout and report with timestamps.  good for
            models that take long time to run.  not good for models that run
            really fast
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
            if platform.system() in 'Windows':
                if not exe_name.lower().endswith('.exe'):
                    exe = which(exe_name + '.exe')
        if exe is None:
            s = 'The program {} does not exist or is not executable.'.format(
                exe_name)
            raise Exception(s)
        else:
            if not silent:
                s = 'pyGSFLOW is using the following executable to run the model: {}'.format(
                    exe)
                print(s)

        exe = os.path.normpath(os.path.join(os.getcwd(), exe))

        if not os.path.isfile(os.path.join(model_ws, namefile)):
            s = 'The namefile for this model does not exists: {}'.format(namefile)
            raise Exception(s)

        # simple little function for the thread to target
        def q_output(output, q):
            for line in iter(output.readline, b''):
                q.put(line)
                # time.sleep(1)
                # output.close()

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
        proc = sp.Popen(argv,
                        stdout=sp.PIPE, stderr=sp.STDOUT, cwd=model_ws)

        if not async_run:
            while True:
                line = proc.stdout.readline()
                c = line.decode('utf-8')
                if c != '':
                    for msg in normal_msg:
                        if msg in c.lower():
                            success = True
                            break
                    c = c.rstrip('\r\n')
                    if not silent:
                        print('{}'.format(c))
                    if report == True:
                        buff.append(c)
                else:
                    break
            return success, buff

        # some tricks for the async stdout reading
        q = Queue.Queue()
        thread = threading.Thread(target=q_output, args=(proc.stdout, q))
        thread.daemon = True
        thread.start()

        failed_words = ["fail", "error"]
        last = datetime.now()
        lastsec = 0.
        while True:
            try:
                line = q.get_nowait()
            except Queue.Empty:
                pass
            else:
                if line == '':
                    break
                line = line.decode().lower().strip()
                if line != '':
                    now = datetime.now()
                    dt = now - last
                    tsecs = dt.total_seconds() - lastsec
                    line = "(elapsed:{0})-->{1}".format(tsecs, line)
                    lastsec = tsecs + lastsec
                    buff.append(line)
                    if not silent:
                        print(line)
                    for fword in failed_words:
                        if fword in line:
                            success = False
                            break
            if proc.poll() is not None:
                break
        proc.wait()
        thread.join(timeout=1)
        buff.extend(proc.stdout.readlines())
        proc.stdout.close()

        for line in buff:
            if normal_msg in line:
                print("success")
                success = True
                break

        if pause:
            input('Press Enter to continue...')
        return success, buff

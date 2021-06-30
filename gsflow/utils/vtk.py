from __future__ import print_function, division
import flopy
import gsflow
import os
import numpy as np

# from flopy.discretization import StructuredGrid


def start_tag(f, tag, indent_level, indent_char="  "):
    s = indent_level * indent_char + tag
    indent_level += 1
    f.write(s + "\n")
    return indent_level


def end_tag(f, tag, indent_level, indent_char="  "):
    indent_level -= 1
    s = indent_level * indent_char + tag
    f.write(s + "\n")
    return indent_level


class Mfvtk(object):
    """
    Generate vtk files for modflow input and output
    """

    def __init__(
        self,
        mf=None,
        vtkname="mf_",
        out_folder=None,
        mf_pkg=[],
        all=True,
        shared_vertex=True,
        ibound_filter=True,
    ):
        if out_folder:
            vtk_nm = os.path.join(out_folder, vtkname)
        else:
            vtk_nm = os.path.join(os.getcwd(), vtkname)
        self.vtkname = vtk_nm
        self.mf = mf
        self.all = all
        if self.all:
            self.mf_pkg = mf.get_package_list()
        else:
            self.mf_pkg = mf_pkg

        self.par3D = ["UPW"]
        self.par2D = ["UZF"]
        self.par1D = ["SFR", "HFB"]
        self.parPoints = ["WEL", "HOB"]

        self.vtk_objects = {}
        self.shared_vertex = shared_vertex
        self.ibound_filter = ibound_filter

    def generate_2d_vtk(self):
        """
        generate vtk objects for each mf package
        :return:
        """
        mf2d = self.__get_dummy_2d_model()
        for pkg in self.mf_pkg:
            if pkg in self.par2D:
                pkg_ = self.mf.get_package(pkg)
                for dataset in pkg_.data_list:
                    if hasattr(dataset, "array"):
                        arr = dataset.array
                        if not (hasattr(arr, "shape")):
                            continue
                        shp = (mf2d.nrow, mf2d.ncol)
                        if arr.shape == shp:
                            if not (pkg in self.vtk_objects.keys()):
                                file_name = self.vtkname + "_" + pkg + ".vtu"
                                obj = Vtk(file_name, mf2d)
                                self.vtk_objects[pkg] = obj
                            arr3d = np.zeros((1, mf2d.nrow, mf2d.ncol))
                            arr3d[0, :, :] = arr
                            nm = dataset.name
                            self.vtk_objects[pkg].add_array(nm, arr3d)

    def generate_1d_vtk(self):
        """
        generate vtk objects for each mf package
        :return:
        """

        for pkg in self.mf_pkg:
            if pkg in self.par1D:
                pkg_ = self.mf.get_package(pkg)
                ibound3d = np.zeros_like(self.mf.bas6.ibound.array)
                rows = pkg_.reach_data["i"]
                cols = pkg_.reach_data["j"]
                lays = pkg_.reach_data["k"]
                ibound3d[lays, rows, cols] = 1
                mf3d = self.__get_dummy_3d_model(ibound=ibound3d)

                columns = pkg_.reach_data.dtype.names
                for field in columns:
                    arr3d = np.zeros((mf3d.nlay, mf3d.nrow, mf3d.ncol))
                    arr3d[lays, rows, cols] = pkg_.reach_data[field]
                    if not (pkg in self.vtk_objects.keys()):
                        file_name = self.vtkname + "_" + pkg + ".vtu"
                        obj = Vtk(file_name, mf3d)
                        self.vtk_objects[pkg] = obj

                    nm = field
                    self.vtk_objects[pkg].add_array(nm, arr3d)

    def __get_dummy_2d_model(self, ibound=None):
        """
        generate 2 dimensional flopy object to be used in
        ploting 2d data
        :return:
        """
        dis = self.mf.dis
        mf2 = flopy.modflow.Modflow("xx")
        nrow = dis.nrow
        ncol = dis.ncol
        delc = dis.delc
        delr = dis.delr
        top = dis.top.array
        botm = top - 0.01
        dis2 = flopy.modflow.ModflowDis(
            mf2,
            nlay=1,
            nrow=nrow,
            ncol=ncol,
            nper=1,
            delr=delr,
            delc=delc,
            top=top,
            botm=botm,
        )
        if not (ibound is None):
            ib = ibound
        else:
            ib = self.mf.bas6.ibound[0, :, :]
        bas = flopy.modflow.ModflowBas(mf2, ibound=ib, strt=1)

        return mf2

    def __get_dummy_3d_model(self, ibound=None):
        """
        generate 2 dimensional flopy object to be used in
        ploting 3d data
        :return:
        """
        dis = self.mf.dis
        mf3 = flopy.modflow.Modflow("xx")
        nrow = dis.nrow
        ncol = dis.ncol
        nlay = dis.nlay
        delc = dis.delc
        delr = dis.delr
        top = dis.top.array
        botm = dis.botm.array
        dis3 = flopy.modflow.ModflowDis(
            mf3,
            nlay=nlay,
            nrow=nrow,
            ncol=ncol,
            nper=1,
            delr=delr,
            delc=delc,
            top=top,
            botm=botm,
        )
        if not (ibound is None):
            ib = ibound
        else:
            ib = self.mf.bas6.ibound.array
        bas = flopy.modflow.ModflowBas(mf3, ibound=ib, strt=1)

        return mf3

    def generate_3d_vtk(self):
        """
        generate vtk objects for each mf package
        :return:
        """
        for pkg in self.mf_pkg:
            if pkg in self.par3D:
                pkg_ = self.mf.get_package(pkg)
                for dataset in pkg_.data_list:
                    if hasattr(dataset, "array"):
                        arr = dataset.array
                        if not (hasattr(arr, "shape")):
                            continue
                        if arr.shape == self.mf.modelgrid.shape:
                            if not (pkg in self.vtk_objects.keys()):
                                file_name = self.vtkname + "_" + pkg + ".vtu"
                                obj = Vtk(file_name, self.mf)
                                self.vtk_objects[pkg] = obj
                            nm = dataset.name[0]
                            self.vtk_objects[pkg].add_array(nm, arr)

    @staticmethod
    def mf_to_vtk(
        vtkname="mf_",
        mfname=None,
        mf_pkg=[],
        all=True,
        out_folder=None,
        shared_vertex=True,
        ibound_filter=True,
    ):

        """

        :param name: modflow name file
        :param packages: a list of packages to visulize
        :param vtkfname: vtk output base file name
        :param all: convert all model's inputs and outputs
        :return:
        """

        # upload the model
        load_only = None
        if len(mf_pkg) > 0:
            load_only = ["DIS", "BAS6"]
            for pkg in mf_pkg:
                if pkg in ["cbc", "hds", "sfr.out"]:
                    continue
                if pkg in load_only:
                    continue
                load_only.append(pkg)

        if not (os.path.isfile(mfname)):
            raise ValueError("{} does not exist...!".format(str(mfname)))

        mf = flopy.modflow.Modflow.load(mfname, load_only=load_only)

        # Write vtk for each package
        if out_folder:
            vtk_nm = os.path.join(out_folder, vtkname)
        else:
            vtk_nm = os.path.join(os.getcwd(), vtkname)

        mfvtk = Mfvtk(
            mf=mf,
            vtkname=vtk_nm,
            mf_pkg=[],
            all=all,
            shared_vertex=shared_vertex,
            ibound_filter=ibound_filter,
        )
        mfvtk.write_vtk()

    def write_vtk(self):

        self.generate_3d_vtk()
        self.generate_2d_vtk()
        self.generate_1d_vtk()
        # mfvtk.generate_0d_vtk()
        self.__write()

    def __write(self):

        for obj in self.vtk_objects.keys():
            shared_vertex = self.shared_vertex
            ibound_filter = self.ibound_filter
            self.vtk_objects[obj].write(
                shared_vertex=shared_vertex, ibound_filter=ibound_filter
            )


class Gsflowvtk:
    def __init__(
        self,
        gs=None,
        vtkname="gs_",
        mf_pkg=[],
        mfall=True,
        out_folder=None,
        shared_vertex=True,
        ibound_filter=True,
    ):
        if out_folder:
            vtk_nm = os.path.join(out_folder, vtkname)
        else:
            vtk_nm = os.path.join(os.getcwd(), vtkname)

        self.vtkname = vtk_nm
        self.gs = gs
        self.mfall = mfall
        self.mf_pkg = mf_pkg
        self.out_folder = out_folder

        self.vtk_objects = {}
        self.shared_vertex = shared_vertex
        self.ibound_filter = ibound_filter

    @staticmethod
    def gsflow_to_vtk(
        vtkname="gs_",
        control_file=None,
        mf_pkg=[],
        mfall=True,
        out_folder=None,
        shared_vertex=True,
        ibound_filter=True,
    ):

        if len(mf_pkg) == 0:
            load_only = None
        else:
            load_only = mf_pkg

        gs = gsflow.GsflowModel.load_from_file(
            control_file, mf_load_only=load_only
        )

        gsfv = Gsflowvtk(
            gs=gs,
            vtkname=vtkname,
            mf_pkg=mf_pkg,
            mfall=mfall,
            out_folder=out_folder,
            shared_vertex=shared_vertex,
            ibound_filter=ibound_filter,
        )
        gsfv.write_vtk()

    def write_vtk(self):
        mf = self.gs.mf
        prms = self.gs.prms

        if mf:
            mfv = Mfvtk(
                mf=mf,
                vtkname=self.vtkname,
                mf_pkg=self.mf_pkg,
                all=self.mfall,
                out_folder=self.out_folder,
                shared_vertex=self.shared_vertex,
                ibound_filter=self.ibound_filter,
            )
            mfv.write_vtk()

        if prms:
            self.write_prms_vtk()

    def write_prms_vtk(self):
        # dummy 2d modflow
        gs = self.gs
        dis = gs.mf.dis
        mf2 = flopy.modflow.Modflow("xx")
        nrow = dis.nrow
        ncol = dis.ncol
        delc = dis.delc
        delr = dis.delr
        top = dis.top.array
        botm = top - 0.01

        dis2 = flopy.modflow.ModflowDis(
            mf2,
            nlay=1,
            nrow=nrow,
            ncol=ncol,
            nper=1,
            delr=delr,
            delc=delc,
            top=top,
            botm=botm,
        )
        hru_type = gs.prms.parameters.get_record("hru_type")
        # TODO:change ibound to hru type

        bas = flopy.modflow.ModflowBas(
            mf2, ibound=hru_type.values.reshape(nrow, ncol), strt=1
        )
        nm = self.vtkname + "_prms.vtu"
        vtkfile = Vtk(nm, mf2)

        # add 2d data
        all_prms = dict()
        gsflow = gs
        nhru = gsflow.prms.parameters.get_record("nhru").values[0]
        for param in gsflow.prms.parameters.record_names:
            par = gsflow.prms.parameters.get_record(param)
            if par.section == "Dimensions":
                continue
            if "nhru" in par.dimensions_names:
                if par.ndim > 1:
                    other_dim = list(np.copy(par.dims))
                    other_dim.remove(nhru)
                    parvalues = par.values.reshape((other_dim[0], nhru))
                    for ipar in range(other_dim[0]):
                        ivar = parvalues[ipar, :]
                        ivar = ivar.reshape(nrow, ncol)  # np.flipud()
                        var = np.zeros((1, nrow, ncol))
                        var[0, :, :] = ivar

                        nm = par.name + "_" + str(ipar + 1)
                        all_prms[nm] = var
                        vtkfile.add_array(param, var)
                else:
                    parvalues = par.values.reshape(nrow, ncol)
                    var = np.zeros((1, nrow, ncol))
                    var[0, :, :] = parvalues

                    all_prms[param] = var
                    vtkfile.add_array(param, var)

        vtkfile.write(
            shared_vertex=self.shared_vertex, ibound_filter=self.ibound_filter
        )


class Vtk(object):
    """
    Support for writing a model to a vtk file
    """

    def __init__(self, output_filename, model, verbose=None):

        if verbose is None:
            verbose = model.verbose
        self.verbose = verbose

        self.output_filename = output_filename
        self.model = model
        self.modelgrid = model.modelgrid
        self.shape = (
            self.modelgrid.nlay,
            self.modelgrid.nrow,
            self.modelgrid.ncol,
        )

        self.arrays = {}

        return

    def add_array(self, name, a):
        assert a.shape == self.shape
        self.arrays[name] = a
        return

    def write(self, shared_vertex=False, ibound_filter=False, htop=None):
        """
        Parameters
        ----------
        shared_vertex : bool
            Make a smoothed representation of model layers by calculating
            an interpolated z value for the cell corner vertices.
        ibound_filter : bool
            Use the ibound array in the basic package of the model that
            was passed in to exclude cells in the vtk file that have an
            ibound value of zero.
        htop : ndarray
            This array must of shape (nlay, nrow, ncol).  If htop is passed
            then these htop values will be used to set the z elevation for the
            cell tops.  This makes it possible to show cells based on the
            saturated thickness.  htop should be calculated by the user as the
            minimum of the cell top and the head and the maximum of the cell
            bottom and the head.
        """
        output_filename = self.output_filename
        assert output_filename.lower().endswith(".vtu")

        if os.path.exists(output_filename):
            if self.verbose:
                print("removing existing vtk file: " + output_filename)
            os.remove(output_filename)

        indent_level = 0
        if self.verbose:
            print("writing vtk file")
        f = open(self.output_filename, "w")

        ibound = None
        if ibound_filter:
            ibound = self.modelgrid.idomain

        dis = self.model.dis
        z = np.vstack(
            [dis.top.array.reshape(1, dis.nrow, dis.ncol), dis.botm.array]
        )
        if shared_vertex:
            verts, iverts = dis.sr.get_3d_shared_vertex_connectivity(
                dis.nlay, z, ibound=ibound
            )
        else:
            top = z[:-1]
            bot = z[1:]
            if htop is not None:
                top = htop
            verts, iverts = dis.sr.get_3d_vertex_connectivity(
                dis.nlay, top, bot, ibound=ibound
            )
        ncells = len(iverts)
        npoints = verts.shape[0]
        if self.verbose:
            s = "Number of point is {}\n " "Number of cells is {}\n".format(
                npoints, ncells
            )
            print(s)

        # xml
        s = '<?xml version="1.0"?>'
        f.write(s + "\n")
        indent_level = start_tag(
            f, '<VTKFile type="UnstructuredGrid">', indent_level
        )

        # unstructured grid
        indent_level = start_tag(f, "<UnstructuredGrid>", indent_level)

        # piece
        s = '<Piece NumberOfPoints="{}" ' 'NumberOfCells="{}">'.format(
            npoints, ncells
        )
        indent_level = start_tag(f, s, indent_level)

        # points
        s = "<Points>"
        indent_level = start_tag(f, s, indent_level)

        s = '<DataArray type="Float64" NumberOfComponents="3">'
        indent_level = start_tag(f, s, indent_level)
        # assert (isinstance(self.modelgrid, StructuredGrid))
        z = np.vstack(
            [
                self.modelgrid.top.reshape(
                    1, self.modelgrid.nrow, self.modelgrid.ncol
                ),
                self.modelgrid.botm,
            ]
        )

        for row in verts:
            s = indent_level * "  " + "{} {} {} \n".format(*row)
            f.write(s)
        s = "</DataArray>"
        indent_level = end_tag(f, s, indent_level)

        s = "</Points>"
        indent_level = end_tag(f, s, indent_level)

        # cells
        s = "<Cells>"
        indent_level = start_tag(f, s, indent_level)

        s = '<DataArray type="Int32" Name="connectivity">'
        indent_level = start_tag(f, s, indent_level)
        for row in iverts:
            s = indent_level * "  " + " ".join([str(i) for i in row]) + "\n"
            f.write(s)
        s = "</DataArray>"
        indent_level = end_tag(f, s, indent_level)

        s = '<DataArray type="Int32" Name="offsets">'
        indent_level = start_tag(f, s, indent_level)
        icount = 0
        for row in iverts:
            icount += len(row)
            s = indent_level * "  " + "{} \n".format(icount)
            f.write(s)
        s = "</DataArray>"
        indent_level = end_tag(f, s, indent_level)

        s = '<DataArray type="UInt8" Name="types">'
        indent_level = start_tag(f, s, indent_level)
        for row in iverts:
            s = indent_level * "  " + "{} \n".format(11)
            f.write(s)
        s = "</DataArray>"
        indent_level = end_tag(f, s, indent_level)

        s = "</Cells>"
        indent_level = end_tag(f, s, indent_level)

        # add cell data
        s = '<CellData Scalars="scalars">'
        indent_level = start_tag(f, s, indent_level)

        self._write_data_array(f, indent_level, "top", z[0:-1], ibound)

        for name, a in self.arrays.items():
            self._write_data_array(f, indent_level, name, a, ibound)

        s = "</CellData>"
        indent_level = end_tag(f, s, indent_level)

        # end piece
        indent_level = end_tag(f, "</Piece>", indent_level)

        # end unstructured grid
        indent_level = end_tag(f, "</UnstructuredGrid>", indent_level)

        # end xml
        indent_level = end_tag(f, "</VTKFile>", indent_level)

        # end file
        f.close()
        return

    def _write_data_array(self, f, indent_level, name, a, ibound):
        """
        Write a numpy array to the vtk file
        """

        # header tag
        s = '<DataArray type="Float64" Name="{}" format="ascii">'.format(name)
        indent_level = start_tag(f, s, indent_level)

        # data
        nlay = a.shape[0]

        # combine ibound with laycbd when model supports laycbd
        if (
            ibound is not None
            and hasattr(self.model, "dis")
            and hasattr(self.model.dis, "laycbd")
        ):
            cbd = np.where(self.model.dis.laycbd.array > 0)
            ibound = np.insert(
                ibound, cbd[0] + 1, ibound[cbd[0], :, :], axis=0
            )

        for k in range(nlay):
            s = indent_level * "  "
            f.write(s)
            if ibound is None:
                ak = a[k].flatten()
            else:
                idx = ibound[k] != 0
                ak = a[k][idx].flatten()
            for v in ak:
                s = " {}".format(v)
                f.write(s)
            f.write("\n")

        # ending tag
        s = "</DataArray>"
        indent_level = end_tag(f, s, indent_level)
        return

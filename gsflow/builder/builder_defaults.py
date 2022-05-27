try:
    import simplejson as json
except ImportError:
    import json

import os


class _DefaultsBase:
    """
    Base class for all defaults classes. Not to be instantiated directly!
    """

    def __init__(self, f=None):
        if f is None:
            ws = os.path.abspath(os.path.dirname(__file__))
            f = os.path.join(ws, "defaults.json")

        self._f = f
        self._json_dict = {}
        self._dict = {}
        self.load_defaults()

    def load_defaults(self):
        """
        Json load method
        """
        with open(self._f) as foo:
            self._json_dict = json.load(foo)

    def to_dict(self, to_json=False):
        raise NotImplementedError()

    def write_defaults_file(self, f):
        """
        Method to write a new JSON defaults file

        Parameters
        ----------
        f : str
            json file name

        """
        d = self.to_dict(to_json=True)
        with open(f, "w") as foo:
            json.dump(d, foo, sort_keys=True, indent=4 * " ")


class Defaults(_DefaultsBase):
    """
    Defaults class that reads and formats JSON based default files
    into an object oriented interface.

    Parameters
    ----------
    f : None, str
        json file name, if None reads in pyGSFLOW's default JSON parameters

    """

    def __init__(self, f=None):
        super().__init__(f)
        self.modflow = ModflowDefaults(f)
        self.prms = PrmsDefaults(f)
        self.control = ControlFileDefaults(f)

    def _create_records(self):
        self.modflow._create_records()
        self.prms._create_records()
        self.control._create_records()

    @property
    def record_names(self):
        """
        Method to get a list of record names

        """
        l = self.modflow.record_names
        l += self.prms.record_names
        l += self.control.record_names
        return l

    def to_dict(self, to_json=False):
        """
        Method to transform the data into a plain dictionary

        Returns
        -------
            dict
        """
        d = {}
        for k, v in self.modflow.to_dict(to_json).items():
            d[k] = v
        for k, v in self.prms.to_dict(to_json).items():
            d[k] = v
        for k, v in self.control.to_dict(to_json).items():
            d[k] = v
        return d


class PrmsDefaults(_DefaultsBase):
    """
    Defaults class that reads and formats JSON based default files
    into an object oriented interface for PRMS.

    Parameters
    ----------
    f : None, str
       json file name, if None reads in pyGSFLOW's default JSON parameters

    """

    def __init__(self, f):
        super().__init__(f)
        self._create_records()

    def _create_records(self):
        if "parameter" in self._json_dict:
            prms = self._json_dict["parameter"]
        else:
            return

        if "dimensions" in prms:
            self._dict["dimensions"] = {}
            for name, value in prms["dimensions"].items():
                self._dict["dimensions"][name] = _DefaultRecord(name, value)

        if "parameters" in prms:
            self._dict["parameters"] = {}
            for name, d in prms["parameters"].items():
                self._dict["parameters"][name] = _PrmsDefaultRecord(
                    name, d["dtype"], d["record"], d["dimension"]
                )

    @property
    def record_names(self):
        """
        Method to get a list of record names

        """
        rnames = []
        for _, d in self._dict.items():
            for name in d.keys():
                rnames.append(name)
        return rnames

    def get_default(self, name):
        """
        Method to get a default parameter object

        Parameters
        ----------
        name: str

        Returns
        -------
        _PrmsDefaultRecord object

        """
        if "dimensions" in self._dict:
            if name in self._dict["dimensions"]:
                return self._dict["dimensions"][name]

        if "parameters" in self._dict:
            if name in self._dict["parameters"]:
                return self._dict["parameters"][name]

    def add_default(self, name, data, dtype=None, dimension=None):
        """
        Method to add or edit an existing default record

        Parameters
        ----------
        name : str
            parameter/dimension name
        data : list, int, or np.ndarray
            list, int, or numpy array of data
        dtype : int or None
            prms dtype, None if it's in the dimension block
        dimension : str or None
            prms dimension type. ex. nssr, nhru, one, etc...

        """
        if dimension is None:
            # assume that this is part of the dimensions block
            rec = _DefaultRecord(name, data)
            if "dimensions" in self._dict:
                self._dict["dimensions"][name] = rec
            else:
                self._dict["dimensions"] = {name: rec}

        else:
            rec = _PrmsDefaultRecord(name, dtype, data, dimension)
            if "parameters" in self._dict:
                self._dict["parameters"][name] = rec
            else:
                self._dict["parameters"] = {name: rec}

    def delete_default(self, name):
        """
        Method to remove a default from the Defaults object

        Parameters
        ----------
        name : str
            parameters/dimension name

        """
        if "dimensions" in self._dict:
            if name in self._dict["dimensions"]:
                self._dict["dimensions"].pop(name)

        if "parameters" in self._dict:
            if name in self._dict["parameters"]:
                return self._dict["parameters"].pop(name)

    def to_dict(self, to_json=False):
        """
        Method to transform the data into a plain dictionary

        Returns
        -------
            dict
        """
        d = {}
        for k, dd in self._dict.items():
            d[k] = {}
            for name, v in dd.items():
                d[k][name] = v.record(to_json)[-1]
        return {"parameter": d}


class ControlFileDefaults(_DefaultsBase):
    """
    Defaults class that reads and formats JSON based default files
    into an object oriented interface for GSFLOW control file.

    Parameters
    ----------
    f : None, str
       json file name, if None reads in pyGSFLOW's default JSON parameters

    """

    def __init__(self, f=None):
        super().__init__(f)
        self._create_records()

    def _create_records(self):
        if "control" in self._json_dict:
            control = self._json_dict["control"]
        else:
            return

        for name, d in control.items():
            self._dict[name] = _PrmsDefaultRecord(
                name, d["dtype"], d["record"]
            )

    @property
    def record_names(self):
        """
        Method to get a list of record names

        """
        return list(self._dict.keys())

    def get_default(self, name):
        """
        Method to get a default parameter object

        Parameters
        ----------
        name: str

        Returns
        -------
        _DefaultRecord object

        """
        if name in self._dict:
            return self._dict[name]

    def add_default(self, name, dtype, record):
        """
        Method to add or edit an existing default record

        Parameters
        ----------
        name : str
            parameter/dimension name
        dtype : int
            gsflow dtype
        record : int, float, str, list, np.ndarray
        """
        rec = _PrmsDefaultRecord(name, dtype, record)
        self._dict[name] = rec

    def delete_default(self, name):
        """
        Method to remove a default from the Defaults object

        Parameters
        ----------
        name : str
            default record name

        """
        if name in self._dict:
            self._dict.pop(name)

    def to_dict(self, to_json=False):
        """
        Method to transform the data into a plain dictionary

        Returns
        -------
            dict
        """
        d = {}
        for name, v in self._dict.items():
            d[name] = v.record(to_json)[-1]
        return {"control": d}


class ModflowDefaults(_DefaultsBase):
    """
    Defaults class that reads and formats JSON based default files
    into an object oriented interface for FloPy modflow.

    Parameters
    ----------
    f : None, str
       json file name, if None reads in pyGSFLOW's default JSON parameters

    """

    def __init__(self, f):
        super().__init__(f)
        self._create_records()

    def _create_records(self):
        pkg_list = []
        for k in self._json_dict.keys():
            if k not in ("control", "parameter"):
                pkg_list.append(k)

        for pkg_name in pkg_list:
            self._dict[pkg_name] = {}
            pkg = self._json_dict[pkg_name]
            for name, value in pkg.items():
                # check for sfr
                if pkg_name.lower() == "sfr":
                    self._dict[pkg_name][name] = {}
                    for k, v in value.items():
                        self._dict[pkg_name][name][k] = _DefaultRecord(k, v)

                else:
                    self._dict[pkg_name][name] = _DefaultRecord(name, value)

    @property
    def record_names(self):
        """
        Method to get a list of record names

        """
        rnames = []
        for pkg_name, pkg in self._dict.items():
            for k, v in pkg.items():
                if not isinstance(v, _DefaultRecord):
                    for kk, vv in v.items():
                        rnames.append(kk)
                else:
                    rnames.append(k)
        return rnames

    def get_default(self, pkg_name, name):
        """
        Method to get a default parameter object

        Parameters
        ----------
        name: str

        Returns
        -------
        _PrmsDefaultRecord object

        """
        pkg_name = pkg_name.lower()
        name = name.lower()
        if pkg_name in self._dict:
            if pkg_name != "sfr":
                if name in self._dict[pkg_name]:
                    return self._dict[pkg_name][name]
            else:
                if name in self._dict[pkg_name]["pkg"]:
                    return self._dict[pkg_name]["pkg"][name]

                elif name in self._dict[pkg_name]["reach"]:
                    return self._dict[pkg_name]["reach"][name]

                elif name in self._dict[pkg_name]["segment"]:
                    return self._dict[pkg_name]["segment"][name]

    def add_default(self, pkg_name, name, data, sfr_block=None):
        """
        Method to add or edit an existing default record

        Parameters
        ----------
        pkg_name : str
            modflow three leter package name
        name : str
            parameter/dimension name
        data : list, int, or np.ndarray
            list, int, or numpy array of data
        sfr_block : str
            if pkg is sfr use sfr block to set to 'pkg', 'reach', or 'segment'

        """
        rec = _DefaultRecord(name, data)
        pkg_name = pkg_name.lower()
        name = name.lower()
        if pkg_name in self._dict:
            if sfr_block is None:
                self._dict[pkg_name][name] = rec
            else:
                self._dict[pkg_name][sfr_block][name] = rec
        else:
            self._dict[pkg_name] = {name: rec}

    def delete_default(self, pkg_name, name=None):
        """
        Method to remove a default from the Defaults object

        Parameters
        ----------
        name : str
            default record name

        """
        pkg_name = pkg_name.lower()
        name = name.lower()
        if pkg_name in self._dict:
            if name is None:
                self._dict.pop(pkg_name)
            else:
                if pkg_name != "sfr":
                    if name in self._dict[pkg_name]:
                        self._dict[pkg_name].pop(name)
                else:
                    if name in self._dict[pkg_name]["pkg"]:
                        self._dict[pkg_name]["pkg"].pop(name)

                    elif name in self._dict[pkg_name]["reach"]:
                        self._dict[pkg_name]["reach"].pop(name)

                    elif name in self._dict[pkg_name]["segment"]:
                        self._dict[pkg_name]["segment"].pop(name)

    def to_dict(self, to_json=False):
        """
        Method to transform the data into a plain dictionary

        Returns
        -------
            dict
        """
        d = {}
        for pkg_name, pkg in self._dict.items():
            d[pkg_name] = {}
            for k, v in pkg.items():
                if not isinstance(v, _DefaultRecord):
                    d[pkg_name][k] = {}
                    for kk, vv in v.items():
                        d[pkg_name][k][kk] = vv.record(to_json)[-1]
                else:
                    d[pkg_name][k] = v.record(to_json)[-1]
        return d


class _DefaultRecord:
    def __init__(self, name, data):
        if name == "stress_period_data":
            if isinstance(data, dict):
                d = {}
                # need to replace string with tuple, OC package
                for k, v in data.items():
                    if isinstance(k, str):
                        x = k.split(",")
                        kstp = int(x[0][1:])
                        kper = int(x[1][:-1])
                        d[(kstp, kper)] = v
                    else:
                        d[k] = v
                data = d

        self.name = name
        self.data = data

    def record(self, to_json=False):
        """
        Record for writing

        Returns
        -------
        str, value
            name and dict of values associted with the name
        """
        if to_json:
            if self.name == "stress_period_data":
                data = {}
                for k, v in self.data.items():
                    if isinstance(k, (tuple, list)):
                        k = str(k)
                    data[k] = v

                return self.name, data

        return self.name, self.data


class _PrmsDefaultRecord(_DefaultRecord):
    def __init__(self, name, dtype, data, dimension=None):
        super().__init__(name, data)
        self.dtype = dtype
        self.dimension = dimension

    def record(self, to_json=False):
        """
        Record for writing

        Returns
        -------
        str, dict
            name and dict of values associted with the name
        """
        d = {"dtype": self.dtype, "record": self.data}
        if self.dimension is not None:
            d["dimension"] = self.dimension

        return self.name, d

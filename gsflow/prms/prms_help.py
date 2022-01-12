import os


class Helper(object):
    """
    Basic Help class
    The function has two goals:
    1) allow the user to get help on a parameter or a variable
    2) allow pyGSFLOW to run error checkers on parameters values and sizes

    """

    def __init__(self):
        self.__prms_parameter_names = None
        self.__prms_dimension_names = None
        self.__prms_output_variables = None
        self._read_param_doc()
        self._read_var_doc()

    @property
    def prms_parameter_names(self):
        return self.__prms_parameter_names

    @property
    def prms_dimension_names(self):
        return self.__prms_dimension_names

    @property
    def prms_output_variables(self):
        return self.__prms_output_variables

    def _read_param_doc(self):
        """
        Reads the PRMS parameter documentation
        """

        fn = os.path.join(
            os.path.dirname(__file__), r"gsflow_prms.control.par_name"
        )
        with open(fn, "r") as fid:
            # content = fid.readlines()
            # fid.close()
            is_dim_section = False
            is_par_section = False
            # content = iter(content)
            dimensions = {}
            parameters = {}
            while True:
                try:
                    line = fid.readline()
                    if "Bounded" in line:
                        curr_par = dimensions[name.strip()]
                        key, value = line.split(":")
                        curr_par[key] = value
                        dimensions[name.strip()] = curr_par
                        continue
                except:
                    # End of file
                    break

                if (
                    line.strip()
                    == "--------------- DIMENSIONS ---------------"
                ):
                    is_dim_section = True
                    is_par_section = False
                    continue

                if (
                    line.strip()
                    == "--------------- PARAMETERS ---------------"
                ):
                    is_par_section = True
                    is_dim_section = False
                    continue

                if not is_dim_section and not is_par_section:
                    continue

                if is_dim_section:
                    if line.strip() == "":
                        # read three lines
                        line = fid.readline()
                        key, name = line.split(":")  # dim name

                        line = fid.readline()
                        key, value = line.split(":")  # dim value

                        line = fid.readline()
                        key, desc = line.split(":")  # dim Desc

                        dimensions[name.strip()] = {
                            "Value": int(value.strip()),
                            "Desc": desc.strip(),
                        }

                if is_par_section:
                    if line.strip() == "":
                        # read three lines
                        line = fid.readline()
                        key, name = line.split(":")  # dim name
                        curr_par = {}
                        for i in range(12):
                            line = fid.readline()
                            try:
                                key, value = line.split(":")  # dim value
                            except:
                                pass
                            value = value.strip()
                            if key.strip() in ["Ndimen", "Size", "Width"]:
                                value = int(value)
                            if key.strip() == "Dimensions":
                                values = value.split(",")
                                value = []
                                for v in values:
                                    dimname, val = v.split("-")
                                    value.append((dimname, int(val)))

                            if key.strip() in ["Max", "Min", "Default"]:
                                if curr_par["Type"] == "float":
                                    value = float(value)
                                elif curr_par["Type"] == "long":
                                    value = int(value)
                                else:
                                    pass  # unknow type

                            curr_par[key.strip()] = value

                        parameters[name.strip()] = curr_par

        self.__prms_parameter_names = parameters
        self.__prms_dimension_names = dimensions

    def _read_var_doc(self):
        """
        Reads the PRMS output variables documentation
        """
        fn = os.path.join(
            os.path.dirname(__file__), r"gsflow_prms.control.var_name"
        )
        fid = open(fn, "r")
        content = fid.read()
        fid.close()
        content = content.split("****")[1]
        content = content[2:]
        content = content.split("\n\n")
        All_variables = dict()
        for line in content:
            recs = line.split("\n")
            curr_record = dict()
            for rec in recs:
                vv = rec.strip().split(":")
                name = vv[0]
                value = " ".join(vv[1:])
                if name in ["Ndimen", "Size"]:
                    value = int(value)

                if name in ["Dimensions"]:
                    if curr_record["Ndimen"] > 1:
                        values = value.split(",")
                        value = []
                        for val in values:
                            val = val.split("-")
                            val[1] = int(val[1])
                            value.append(val)
                    else:
                        value = value.split("-")
                        value[1] = int(value[1])

                if name == "Name":
                    Field_name = value
                else:
                    curr_record[name] = value
            All_variables[Field_name] = curr_record

            self.__prms_output_variables = All_variables

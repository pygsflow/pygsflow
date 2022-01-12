class GsConstant(object):
    """
    Class object that holds global constants
    """

    COLUMN_HEADER = ["Year", "Month", "Day", "Hour", "Minute", "Second"]
    PRMS_DATA_TYPES = {
        1: "integer",
        2: "real",
        3: "double",
        4: "string",
        "4": "string",
    }
    PRMS_SECTIONS = ("Dimensions", "Parameters")
    GSFLOW_FILES = [
        "csv_output_file",
        "data_file",
        "gsflow_output_file",
        "model_output_file",
        "modflow_name",
        "param_file",
        "stat_var_file",
        "var_init_file",
        "var_save_file",
        "ani_output_file",
        "param_print_file",
        "stats_output_file",
        "nsubOutBaseFileName",
        "basinOutBaseFileName",
        "nsegmentOutBaseFileName",
        "nhruOutBaseFileName",
    ]

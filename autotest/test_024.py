import os
import numpy as np

ws = os.path.abspath(os.path.dirname(__file__))
cellsize = 90


def test_import_prms_builder():
    from gsflow.builder import PrmsBuilder


def test_build_prms_parameters():
    from gsflow.builder import GenerateFishnet, FlowAccumulation, PrmsBuilder
    from gsflow.prms import PrmsParameters

    param_out = os.path.join(ws, "temp", "builder.prms.parameter")
    data_ws = os.path.join(ws, "..", "examples", "data", "geospatial")
    resampled_ws = os.path.join(data_ws, "resampled_90m")
    dem_file = os.path.join(data_ws, "dem.img")
    resampled_dem = os.path.join(resampled_ws, "dem_90m_median.txt")
    watershed = os.path.join(resampled_ws, "watershed_90m.txt")
    streams = os.path.join(resampled_ws, "streams_90m.bin")
    cascades = os.path.join(resampled_ws, "cascades_90m.bin")

    modelgrid = GenerateFishnet(
        dem_file,
        xcellsize=cellsize,
        ycellsize=cellsize
    )

    dem_data = np.genfromtxt(resampled_dem)
    watershed = np.genfromtxt(watershed, dtype=int)
    streams = FlowAccumulation.load_streams(streams)
    cascades = FlowAccumulation.load_cascades(cascades)

    prms_build = PrmsBuilder(
        streams,
        cascades,
        modelgrid,
        dem_data,
        watershed
    )

    parameters = prms_build.build("sagehen")

    if not isinstance(parameters, PrmsParameters):
        raise AssertionError("PrmsBuilder returned wrong type")

    parameters.write(param_out)

    parameters2 = PrmsParameters.load_from_file(param_out)

    for record in parameters.record_names:
        valid = parameters.get_values(record)
        test = parameters.get_values(record)

        if not np.allclose(valid, test):
            raise AssertionError("PrmsBuilder create or write error")


if __name__ == "__main__":
    test_import_prms_builder()
    test_build_prms_parameters()
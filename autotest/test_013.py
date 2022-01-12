# plotting testing
import os
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.tri.tricontour import TriContourSet
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from gsflow import GsflowModel
from gsflow.output import PrmsDiscretization, PrmsPlot


ws = os.path.abspath(os.path.dirname(__file__))


def test_prms_plotting():
    local_ws = os.path.join(ws, "..", "examples", "data", "sagehen", "gsflow")
    control_file = "saghen_new_cont.control"
    gs = GsflowModel.load_from_file(os.path.join(local_ws, control_file))

    shp_ws = os.path.join(ws, "..", "examples", "data", "sagehen", "shapefiles")
    shp = "hru_params.shp"

    dis = PrmsDiscretization.load_from_shapefile(os.path.join(shp_ws, shp))
    assert isinstance(dis, PrmsDiscretization)
    assert dis.nhru == 6468
    assert dis.extent is not None

    plot = PrmsPlot(dis)

    ssr2gw = gs.prms.parameters.get_record("ssr2gw_rate")

    ax = plot.plot_parameter(ssr2gw, masked_values=[0], cmap="viridis")
    assert isinstance(ax, PatchCollection)
    plt.close()

    ax = plot.contour_parameter(ssr2gw, masked_values=[0])
    assert isinstance(ax, TriContourSet)
    plt.close()

    rain_adj = gs.prms.parameters.get_record("rain_adj")

    fig = plot.plot_parameter(rain_adj, masked_values=[0], cmap="jet")
    assert isinstance(fig, Figure)
    plt.close()

    fig = plot.contour_parameter(rain_adj)
    assert isinstance(fig, Figure)
    plt.close()

    ax = plot.plot_model_discretization()
    assert isinstance(ax, Axes)
    plt.close()

    ax = dis.plot_discretization()
    assert isinstance(ax, Axes)
    plt.close()


if __name__ == "__main__":
    test_prms_plotting()

import os

ws = os.path.abspath(os.path.dirname(__file__))


def test_fishnet_import():
    from gsflow.builder import GenerateFishnet


def test_fishnet_generation():
    from gsflow.builder import GenerateFishnet
    from flopy.discretization import StructuredGrid

    raster = os.path.join(
        ws, '..', 'examples', 'data', 'geospatial', 'dem.img'
    )
    modelgrid = GenerateFishnet(raster, xcellsize=50, ycellsize=100, buffer=0)

    if not isinstance(modelgrid, StructuredGrid):
        raise TypeError('GenerateFishnet should return SturcturedGrid object')

    if modelgrid.ncol != 149 or modelgrid.nrow != 69:
        raise AssertionError("DELC and DELR have not been properly calculated")


if __name__ == "__main__":
    test_fishnet_import()
    test_fishnet_generation()

# test the instantiation of empty classes
import gsflow

def test_import_classes():
    from gsflow import ParameterRecord
    from gsflow import ControlRecord
    from gsflow import PrmsParameters
    from gsflow import PrmsData
    from gsflow import PrmsModel
    from gsflow.modflow import Modflow
    from gsflow import ControlFile
    from gsflow import GsflowModel


if __name__ == "__main__":
    test_import_classes()

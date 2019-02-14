from gsflow.modflow import Modflow


def test_empty_modflow():
    modelname = "gsflow_test"
    ml = Modflow(modelname=modelname)

    assert isinstance(ml, Modflow)
    assert ml.name == modelname
    assert ml.version == 'mfnwt'


if __name__ == "__main__":
    test_empty_modflow()

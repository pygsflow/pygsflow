
def test_import_gsflow():
    success = True
    try:
        import gsflow
    except:
        success = False
    assert success


if __name__ == "__main__":
    test_import_gsflow()

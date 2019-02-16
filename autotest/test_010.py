import gsflow

def test_call_help():
    helper = gsflow.prms.Helper()
    assert(len(helper.prms_dimension_names) > 0)
    assert(len(helper.prms_output_variables) > 0)
    assert(len(helper.prms_parameter_names) > 0)


if __name__ == "__main__":
    test_call_help()
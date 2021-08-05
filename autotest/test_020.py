import os


ws = os.path.abspath(os.path.dirname(__file__))


def test_import_defaults():
    from gsflow.builder import Defaults
    from gsflow.builder import ModflowDefaults
    from gsflow.builder import PrmsDefaults
    from gsflow.builder import ControlFileDefaults


def test_defaults():
    from gsflow.builder import Defaults

    outfile = os.path.join(ws, 'temp', 'test_20.json')

    defaults = Defaults()

    if len(defaults.record_names) != 166:
        raise AssertionError("record_names failed")

    control = defaults.control
    r = control.get_default('rpt_days')
    control.add_default(r.name, r.data + 2)

    r2 = control.get_default('rpt_days')
    if r2.data != 9:
        raise AssertionError("control.add_default failed")

    control.delete_default('rpt_days')

    if 'rpt_days' in control.record_names:
        raise AssertionError("control.delete_default failed")

    prms = defaults.prms
    r = prms.get_default("elev_units")
    prms.add_default(r.name, 0, 1, 'one')
    r2 = prms.get_default(r.name)

    if r2.data != 0:
        raise AssertionError('prms.add_default failed')

    prms.delete_default('elev_unit')

    if 'elev_unit' in prms.record_names:
        raise AssertionError("prms.delete_default failed")

    modflow = defaults.modflow
    r = modflow.get_default('uzf', 'ntrail2')
    modflow.add_default('uzf', r.name, r.data + 10)
    r2 = modflow.get_default('uzf', r.name)

    if r2.data != 20:
        raise AssertionError('modflow.add_default failed')

    modflow.delete_default('uzf', 'ntrail2')

    if 'ntrail2' in modflow.record_names:
        raise AssertionError("modflow.delete_default failed")

    defaults.prms.add_default('elev_units', 0, 1, 'one')
    defaults.prms.delete_default('cov_type')
    defaults.write_defaults_file(outfile)

    defaults2 = Defaults(outfile)

    if len(defaults2.record_names) != 163:
        raise AssertionError("write or read custom file function failed")

    r = defaults2.prms.get_default("elev_units")
    if r.data != 0:
        raise AssertionError("write or read custom file function failed")


if __name__ == "__main__":
    test_import_defaults()
    test_defaults()
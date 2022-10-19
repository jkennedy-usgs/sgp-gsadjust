from gsadjust.plots.datum import PlotDatumCompare


def test_PlotDatumCompare(mainprog, obstreesurvey_with_results):
    test_plot = PlotDatumCompare(obstreesurvey_with_results, mainprog)
    assert test_plot.figure is not None
    assert test_plot.survey.name == "2020-05-11"
    assert test_plot.get_data()[0][0] < 0.001
    assert test_plot.get_data()[1][0] == "rg37"

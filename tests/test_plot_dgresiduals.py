from gsadjust.plots.residual import PlotDgResidualHistogram

def test_PlotDgResidualHistogram(mainprog, obstreesurvey_with_results):
    test_plot = PlotDgResidualHistogram(obstreesurvey_with_results, mainprog)
    assert test_plot.figure is not None
    assert test_plot.get_data()[0] - 4.8 < 0.001
    assert len(test_plot.get_data()) == 45
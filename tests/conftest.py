import sys, os

code_path = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, code_path + '..')
sys.path.insert(0, code_path)

import pytest
import gsadjust
from gsadjust.obstree.survey import ObsTreeSurvey
import pickle
from PyQt5 import QtGui, QtCore, QtWidgets
from gsadjust.gui.tabs import drift


@pytest.fixture
def test_channellist_fixture(request):
    fname = request.param  #'channellist_testobj.p'
    with open(fname, "rb") as f:
        cl = pickle.load(f)
    return cl


@pytest.fixture
def list_of_deltas():
    """
    Deltas can't just be loaded from a pickle file, because they contain ObsTreeStation objects (pyqt objects can't be
    pickled). Instead, load deltas from a workspace. It's a bit circular, because we have to assume the deltas were
    created correctly in the first place.
    :return:
    """
    delta_list = []
    fname = './tests/test_workspace1.gsa'
    # app = GSadjust.MainProg()
    # app.workspace_open(fname)
    ots, dms, _ = gsadjust.MainProg.obsTreeModel.load_workspace(fname)
    survey = ots[0]
    # survey = app.obsTreeModel.invisibleRootItem().child(0)
    loop = survey.child(0, 0)
    delta_model = loop.delta_model
    for i in range(delta_model.rowCount()):
        idx = delta_model.index(i, 0)
        delta_list.append(delta_model.data(idx, role=QtCore.Qt.UserRole))
    return delta_list


@pytest.fixture
def test_twostations_fixture(request):
    """
    For testing the creation of normal deltas.
    :param request: pytest attribute, stores filename to load
    :return: tuple (ObsTreeStation1, ObsTreeStation2)
    """
    fname = request.param
    with open(fname, "rb") as f:
        cl = pickle.load(f)
    obstreesurvey = ObsTreeSurvey('test')
    obstreesurvey.populate(cl, name='0')
    loop = obstreesurvey.child(0)
    return (loop.child(0), loop.child(1))


@pytest.fixture
def test_threestations_fixture(request):
    """
    For testing three-point deltas.
    :param request: pytest attribute, stores filename to load
    :return: tuple (ObsTreeStation1, ObsTreeStation2, ObsTreeStation3)
    """
    fname = request.param[0]
    idxs = request.param[1]
    with open(fname, "rb") as f:
        cl = pickle.load(f)
    obstreesurvey = ObsTreeSurvey('test')
    obstreesurvey.populate(cl, name='0')
    loop = obstreesurvey.child(0)
    return (loop.child(idxs[0]), loop.child(idxs[1]), loop.child(idxs[2]))


@pytest.fixture
def channellist():
    fname = './tests/channellist_burris_single.p'
    with open(fname, "rb") as f:
        cl = pickle.load(f)
    return cl


@pytest.fixture
def obstreestation():
    a = gsadjust.obstree.station(channellist(), 'jeff', 1)
    return a


@pytest.fixture
def obstreesurvey():
    filename = './tests/test_BurrisData.txt'
    meter_type = 'Burris'
    survey = ObsTreeSurvey('2020-05-11')
    assert type(survey) is ObsTreeSurvey
    data = gsadjust.GSadjust.MainProg.read_raw_data_file(filename, meter_type)
    survey.populate(data)
    return survey


@pytest.fixture
def obstreesurvey_with_results(obstreesurvey):
    from gsadjust.data import Datum

    obstreeloop = obstreesurvey.child(0)
    data = obstreeloop.checked_stations()
    model = gsadjust.gui.tabs.TabDrift.calc_none_dg(data, 'test')
    obstreeloop.deltas = model
    survey = obstreesurvey
    survey.populate_delta_model()
    d = Datum('rg37')
    survey.datums.append(d)
    survey.run_inversion()
    return survey


@pytest.fixture
def obstreemodel(obstreesurvey):
    model = gsadjust.GSadjust.ObsTreeModel()
    model.appendRow([obstreesurvey, QtGui.QStandardItem('a'), QtGui.QStandardItem('a')])
    assert type(model) is gsadjust.GSadjust.ObsTreeModel
    return model


@pytest.fixture
def obstreemodel_with_results(obstreesurvey_with_results):
    model = gsadjust.GSadjust.ObsTreeModel()
    model.appendRow(
        [obstreesurvey_with_results, QtGui.QStandardItem('a'), QtGui.QStandardItem('a')]
    )
    assert type(model) is gsadjust.GSadjust.ObsTreeModel
    return model


@pytest.fixture
def mainprog():
    fname = './tests/test_workspace10.gsa'
    # app = QtWidgets.QApplication(sys.argv)
    mp = gsadjust.GSadjust.MainProg()
    mp.workspace_clear(confirm=False)
    mp.workspace_open_json(fname)
    return mp

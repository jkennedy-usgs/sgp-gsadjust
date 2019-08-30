import sys, os
code_path = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, code_path + '/../gsadjust')
sys.path.insert(0, code_path)

import pytest
import pyqt_models
import GSadjust
import pickle
from PyQt5 import QtGui, QtCore



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
    fname = './tests/test_workspace1.p'
    # app = GSadjust.MainProg()
    # app.workspace_open(fname)
    ots, dms, _ = GSadjust.MainProg.obsTreeModel.load_workspace_p(fname)
    survey = ots[0]
    # survey = app.obsTreeModel.invisibleRootItem().child(0)
    loop = survey.child(0,0)
    delta_model = loop.delta_model
    for i in range(delta_model.rowCount()):
        idx = delta_model.index(i,0)
        delta_list.append(delta_model.data(idx, role=QtCore.Qt.UserRole))
    jeff = 1
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
    obstreesurvey = pyqt_models.ObsTreeSurvey('test')
    obstreesurvey.populate(cl, name='test')
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
    obstreesurvey = pyqt_models.ObsTreeSurvey('test')
    obstreesurvey.populate(cl, name='test')
    loop = obstreesurvey.child(0)
    return (loop.child(idxs[0]), loop.child(idxs[1]), loop.child(idxs[2]))

@pytest.fixture
def channellist():
    fname = './tests/channellist_burris.p'
    with open(fname, "rb") as f:
        cl = pickle.load(f)
    return cl

@pytest.fixture
def obstreestation():
    a = pyqt_models.ObsTreeStation(channellist(),'jeff',1)
    return a

@pytest.fixture
def obstreesurvey():
    filename = './tests/test_BurrisData.txt'
    meter_type = 'Burris'
    survey = pyqt_models.ObsTreeSurvey('test')
    assert type(survey) is pyqt_models.ObsTreeSurvey
    data = GSadjust.MainProg.read_raw_data_file(filename, meter_type)
    survey.populate(data)
    return survey

@pytest.fixture
def obstreemodel(obstreesurvey):
    model = pyqt_models.ObsTreeModel()
    model.appendRow([obstreesurvey, QtGui.QStandardItem('a'), QtGui.QStandardItem('a')])
    assert type(model) is pyqt_models.ObsTreeModel
    return model

@pytest.fixture
def mainprog():
    fname = './tests/test_workspace1.p'
    app = GSadjust.MainProg()
    app.workspace_open(fname)
    return app
import pytest
import pyqt_models


@pytest.fixture
def channellist():
    import pickle
    fname = 'channellist_testobj.p'
    with open(fname, "rb") as f:
        cl = pickle.load(f)
    return cl

@pytest.fixture
def obstreestation():
    a = pyqt_models.ObsTreeStation(channellist(),'jeff',1)
    return a

@pytest.fixture
def obstreesurvey():
    a = pyqt_models.ObsTreeSurvey('test')
    return a
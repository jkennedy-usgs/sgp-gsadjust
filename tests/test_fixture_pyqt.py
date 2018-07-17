import pytest
import pyqt_models
import GSadjust

FILES = ['channellist_testobj.p',
         'channellist_testobj.p',
         'channellist_testobj.p']

@pytest.fixture(files=FILES)
def channellist(file):
    import pickle
    fname = file#'channellist_testobj.p'
    with open(fname, "rb") as f:
        cl = pickle.load(f)
    return cl

def scintrex_channellist():
    import pickle
    fname = 'channellist_testobj.p'
    with open(fname, "rb") as f:
        cl = pickle.load(f)
    return cl

def cg6_chanellist():
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
    filename = './test_BurrisData.txt'
    meter_type = 'Burris'
    obstreesurvey = pyqt_models.ObsTreeSurvey('test')
    assert type(obstreesurvey) is pyqt_models.ObsTreeSurvey
    data = GSadjust.MainProg.read_raw_data_file(filename, meter_type)
    obstreesurvey.populate(data)
    return obstreesurvey
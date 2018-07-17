import pytest
import pyqt_models
import GSadjust


@pytest.fixture
def test_channellist_fixture(request):
    import pickle
    fname = request.param  #'channellist_testobj.p'
    with open(fname, "rb") as f:
        cl = pickle.load(f)
    return cl


@pytest.fixture
def test_twostations_fixture(request):
    import pickle
    fname = request.param
    with open(fname, "rb") as f:
        cl = pickle.load(f)
    obstreesurvey = pyqt_models.ObsTreeSurvey('test')
    obstreesurvey.populate(cl, name='test')
    loop = obstreesurvey.child(0)
    return (loop.child(0), loop.child(1))


@pytest.fixture
def channellist():
    import pickle
    fname = 'channellist_burris.p'
    with open(fname, "rb") as f:
        cl = pickle.load(f)
    return cl

@pytest.fixture
def obstreestation():
    a = pyqt_models.ObsTreeStation(channellist(),'jeff',1)
    return a

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
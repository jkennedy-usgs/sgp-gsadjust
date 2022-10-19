import sys
import os
import pytest


def test_numpy_inversion(obstreesurvey):
    import numpy as np

    # Test problem from Adj. Computations Spatial Data Analysis, Ghilani and Wolf, Wiley
    # ftp://doc.nit.ac.ir/civil/m.abbaszadeh/Theory%20of%20Errors%20and%20Adjustment/
    # ebooksclub.org__Adjustment_Computations__Spatial_Data_Analysis.pdf
    A = [[1, 0, 0], [-1, 1, 0], [0, -1, 1], [0, 0, -1], [-1, 0, 1], [0, 1, 0]]
    obstreesurvey.adjustment.A = np.array(A)

    P = [
        [1 / 0.006**2, 0, 0, 0, 0, 0],
        [0, 1 / 0.004**2, 0, 0, 0, 0],
        [0, 0, 1 / 0.005**2, 0, 0, 0],
        [0, 0, 0, 1 / 0.003**2, 0, 0],
        [0, 0, 0, 0, 1 / 0.004**2, 0],
        [0, 0, 0, 0, 0, 1 / 0.012**2],
    ]
    obstreesurvey.adjustment.P = np.array(P)

    obs = [448.105, 5.360, -8.523, -444.944, -3.167, 453.477]
    obstreesurvey.adjustment.Obs = obs

    obstreesurvey.adjustment.dof = 3

    obstreesurvey.adjustment.python_lsq_inversion()

    answer = np.array([448.1087, 453.4685, 444.9436])
    answer_sd = 0.6575
    diff = answer - obstreesurvey.adjustment.X

    assert max(abs(diff)) < 0.0001
    assert obstreesurvey.adjustment.SDaposteriori - answer_sd < 0.0001

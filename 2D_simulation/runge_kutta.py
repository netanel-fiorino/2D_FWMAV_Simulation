import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('..')

from parameters import *
from derivatives import dxdt

def runge_kutta(time_sim, current, tmp_left, tmp_right, delta_t):
    state= current
    K1 = delta_t*dxdt(time_sim, state, tmp_left, tmp_right)
    K2 = delta_t*dxdt(time_sim, state + 1/2*K1, tmp_left, tmp_right)
    K3 = delta_t*dxdt(time_sim, state + 1/2*K2, tmp_left, tmp_right)
    K4 = delta_t*dxdt(time_sim, state + 1*K3, tmp_left, tmp_right)
    state += (K1/6 + K2/3 + K3/3 + K4/6)
    return state
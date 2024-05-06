import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('..')

from parameters import *
from WingShape import WingShape
from tail_force import aero_force_tail
from body_force import aero_force_body

def dxdt(t, state, tmp_left, tmp_right):
    dx_dt = [0, 0, 0, 0, 0, 0]          #xdot, ydot, thetadot, xddot, yddot, thetaddot
    x, y, theta, xdot, ydot, thetadot = state
    fx_aero_tail, fy_aero_tail, alpha_tail = aero_force_tail(t, theta, xdot, ydot, thetadot)
    fx_aero_body, fy_aero_body, alpha_body = aero_force_body(t, theta, xdot, ydot, thetadot)
    fx_aero = fx_aero_tail + fx_aero_body
    fy_aero = fy_aero_tail + fy_aero_body
    pin = False     
    if pin == True:
        dx_dt[0] = 0
        dx_dt[1] = 0
        dx_dt[2] = 0
        dx_dt[3] = 0
        dx_dt[4] = 0
        dx_dt[5] = 0
    else:
        dx_dt[0] = xdot
        dx_dt[1] = ydot
        dx_dt[2] = thetadot
        dx_dt[3] = ((-(tmp_left + tmp_right).item(2) + fy_aero)*-np.sin(theta) + ((tmp_left + tmp_right).item(0) + fx_aero)*np.cos(theta))/m
        dx_dt[4] = (-m*g + (-(tmp_left + tmp_right).item(2) + fy_aero)*np.cos(theta) + ((tmp_left + tmp_right).item(0) + fx_aero)*np.sin(theta))/m
        dx_dt[5] = (r_wing*-(tmp_left + tmp_right).item(2) + r_tail*fy_aero_tail + r_body*fy_aero_body + (tmp_left + tmp_right).item(4))/(Jz*1)
    return np.array(dx_dt)
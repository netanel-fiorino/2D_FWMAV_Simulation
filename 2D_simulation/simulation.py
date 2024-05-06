import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('..')

from parameters import *
from WingShape import WingShape
from WingDynamics import WingDynamics
from tail_force import aero_force_tail
from runge_kutta import runge_kutta

def simulation(x0, y0, angle0, vx0, vy0, omega0, stroke_offset, rotation_offset):
    tArr = np.zeros(int(num_periods*t_period/t_step))
    alpha_Arr = np.zeros(int(num_periods*t_period/t_step))
    results = np.zeros((int(num_periods*t_period/t_step), 6))
    wing_force = np.zeros((int(num_periods*t_period/t_step), 2))
    wing_moment = np.zeros((int(num_periods*t_period/t_step), 1))
    wing_angles = np.zeros((int(num_periods*t_period/t_step), 3))
    wing_angle_vel = np.zeros((int(num_periods*t_period/t_step), 3))

    time_sim = 0
    time_flapping = 0
    i = 0

    Wingshape = WingShape(wing_length, AR_wing, 0)
    # wing_left = WingDynamics(side='left', Wingshape=Wingshape, v0=np.array([[v0*np.cos(np.radians(angle0))], [-v0*np.sin(np.radians(angle0))], [0]]), glide=True, phi=0, phi_dot=0, theta=0, theta_dot=0)               #rectangle ellipse
    # wing_right = WingDynamics(side='right', Wingshape=Wingshape, v0=np.array([[v0*np.cos(np.radians(angle0))], [-v0*np.sin(np.radians(angle0))], [0]]), glide=True, phi=0, phi_dot=0, theta=0, theta_dot=0)

    wing_left = WingDynamics(side='left', Wingshape=Wingshape, v0=np.array([[vx0], [0], [-vy0]]), glide=False)               #rectangle ellipse
    wing_right = WingDynamics(side='right', Wingshape=Wingshape, v0=np.array([[vx0], [0], [-vy0]]), glide=False)


    state = [x0, y0, np.radians(angle0), vx0, vy0, omega0]

    wing_left._state[3][0] = (state[3]*np.cos(state[2]) + state[4]*np.sin(state[2]))
    wing_left._state[4][0] = 0
    wing_left._state[5][0] = -(state[4]*np.cos(state[2]) - state[3]*np.sin(state[2]) )#+ r_wing*state[5])

    wing_right._state[3][0] = (state[3]*np.cos(state[2]) + state[4]*np.sin(state[2]))
    wing_right._state[4][0] = 0
    wing_right._state[5][0] = -(state[4]*np.cos(state[2]) - state[3]*np.sin(state[2]) )#+ r_wing*state[5])

    wing_right.update_initial_angle(-2*t_step, stroke_offset, rotation_offset)
    wing_left.update_initial_angle(-2*t_step, stroke_offset, rotation_offset)
    #print(wing_left.V_W_inflow_prev)
    tmp_left = wing_left.force_calc()
    tmp_right = wing_right.force_calc()

    wing_right.update_initial_angle(-t_step, stroke_offset, rotation_offset)
    wing_left.update_initial_angle(-t_step, stroke_offset, rotation_offset)
    #print(wing_left.V_W_inflow_prev)
    tmp_left = wing_left.force_calc()
    tmp_right = wing_right.force_calc()

    wing_right.update_initial_angle(time_flapping, stroke_offset, rotation_offset)
    wing_left.update_initial_angle(time_flapping, stroke_offset, rotation_offset)
    #print(wing_left.V_W_inflow_prev)
    tmp_left = wing_left.force_calc()
    tmp_right = wing_right.force_calc()
    #print(wing_left.V_W_inflow_prev)

    while round(time_sim, 9) < round(num_periods*t_period, 9):
        results[i, :] = state
        if np.sign(state[2]) == 1:
            state[2] = state[2] % (2 * np.pi)
        else:
            state[2] = state[2] % (-2 * np.pi)

        wing_force[i, :] = [(tmp_left + tmp_right).item(0), -(tmp_left + tmp_right).item(2)]
        wing_moment[i, :] = [(tmp_left + tmp_right).item(4)]
        wing_angles[i, :] = [wing_left._state[7][0], wing_left._state[8][0], wing_left._state[9][0]]
        wing_angle_vel[i, :] = [wing_left._state[10][0], wing_left._state[11][0], wing_left._state[12][0]]

        tArr[i] = time_sim
        alpha_Arr[i] = wing_left.alpha_COP

        time_sim = time_sim + t_step
        if wing_left.glide == False:
            time_flapping = time_flapping + t_step
            if round(time_flapping, 8) >= round(t_period, 8):
                time_flapping = 0



        wing_left._state[3][0] = (state[3]*np.cos(state[2]) + state[4]*np.sin(state[2]))
        wing_left._state[4][0] = 0
        wing_left._state[5][0] = -(state[4]*np.cos(state[2]) - state[3]*np.sin(state[2]))# + r_wing*state[5])

        wing_right._state[3][0] = (state[3]*np.cos(state[2]) + state[4]*np.sin(state[2]))
        wing_right._state[4][0] = 0
        wing_right._state[5][0] = -(state[4]*np.cos(state[2]) - state[3]*np.sin(state[2]))# + r_wing*state[5])

        wing_right.update_initial_angle(time_flapping, stroke_offset, rotation_offset)
        wing_left.update_initial_angle(time_flapping, stroke_offset, rotation_offset)

        tmp_left = wing_left.force_calc()
        tmp_right = wing_right.force_calc()
        #print((state[3]*np.cos(state[2]) + state[4]*np.sin(state[2])))
        state = runge_kutta(time_sim, state, tmp_left, tmp_right, t_step)
        i+=1
    
    tail_force = aero_force_tail(tArr, results[:,2], results[:,3], results[:,4], results[:,5])
    return tArr, results, wing_force, tail_force, alpha_Arr, wing_angles, wing_angle_vel
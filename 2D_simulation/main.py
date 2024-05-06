import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('..')

from parameters import *
from simulation import simulation
from plots import plot_figures

x0 = 0
y0 = 10
angle0 = 0      #initial pitch in degrees
vx0 = 2
vy0 = 0
omega0 = 0      #initial pitch rate in rad/s
stroke_offset = 0       #initial stroke angle
rotation_offset = 0        #inital rotation angle


tArr, results, wing_force, tail_force, alpha_Arr, wing_angles, wing_angle_vel = simulation(x0, y0, angle0, vx0, vy0, omega0, stroke_offset, rotation_offset)

plot_figures(tArr, results, wing_force, tail_force, alpha_Arr, wing_angles, wing_angle_vel)
import numpy as np
import matplotlib.pyplot as plt

#####################################################
                #Physical Parameters
#####################################################
r_wing = 1/3*0.0381          #location from COM to Wing Pivot
r_tail = -0.1            #location from COM to Tail Center of Pressure
m = 0.03                   #Mass [kg]
g = 9.81*1.0                #gravitational constant [m/s^2]
rho = 1.2682                #Density of air [kg/m^3]
L = 0.254                      #Length of FWMAV
H = 0.02                    #Height of FWMAV

angle_tail = np.radians(20)
S_tail = 0.2*0.0677       #surface area of tail
l_tail = 0.1778               #tail length
AR_tail = l_tail**2/S_tail

S_body = 0.00677       #surface area of tail
l_body = 0.05               #tail length
AR_body = l_body**2/S_body
r_body = -0.03

wing_length = 0.152        #m
AR_wing = 3.25

e = 0.9                     #oswald efficiency
M = 25                      #for blending function
alpha0 = 0.0                #for blending function

Jx = 1/12*m*(L**2+H**2)
Jy = 1/12*m*(H**2)
Jz = 0.9* 1/12*m*(L**2)

#####################################################
                #Simulation Parameters
#####################################################
num_elements = 5            #Number of Blade elements per wing

t_period = .05              #sec
freq = 10                   #1/t_period   #Hz
t_period = 1/freq           #sec (if you want to define the speed by frequency)

num_seconds = 2                    #simulation time
num_periods = 10 #num_seconds/t_period

t_step = t_period/200

#####################################################
                #Initial Conditions
#####################################################
v_free_stream = 0*np.array([[1], [0], [0]])     #m/s
beta = np.radians(90)                           #stroke plane angle

phi_start = np.radians(-25)                           #initial stroke angle
phi0 = np.radians(10)                           #Average Stroke angle
phi_amplitude = np.radians(70/2)                   #Amplitude of stroke angle

psi_start = np.radians(0)                               #initial deviation angle/elevation angle
psi0 = np.radians(0)                                    #Average deviation angle
psi_amplitude = np.radians(0)

theta_start = np.radians(-5)                        #initial pitch angle/rotation angle
theta0 = np.radians(-5)                          #Average Pitch angle
theta_amplitude = np.radians(15/2)                        #amplitude of pitch rotation

C_alpha = 2.6
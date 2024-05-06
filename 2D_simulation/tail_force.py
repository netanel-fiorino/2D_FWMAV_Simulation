import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('..')

from parameters import *
from WingShape import WingShape

def aero_force_tail(t, theta, xdot, ydot, thetadot):
    # Ensure that all inputs are NumPy arrays
    t = np.array(t)
    theta = np.array(theta)
    xdot = np.array(xdot)
    ydot = np.array(ydot)
    thetadot = np.array(thetadot)

    # Calculate tail_xdot and tail_ydot for arrays
    tail_ydot = ydot + thetadot * r_tail * np.cos(theta)
    tail_xdot = xdot + thetadot * r_tail * -np.sin(theta)
    #print(f"velocity: {tail_xdot} and {tail_ydot}")
    # Handle the case when theta is an array
    if theta.size > 1:
        theta = np.where(theta > 0, theta % (2 * np.pi), theta % (-2 * np.pi))
    else:
        if np.sign(theta) == 1:
            theta = theta % (2 * np.pi)
        else:
            theta = theta % (-2 * np.pi)

    alpha = theta - np.arctan2(tail_ydot, tail_xdot) - angle_tail
    if alpha.size > 1:
        alpha = np.where(alpha > 0, alpha % (2 * np.pi), alpha % (-2 * np.pi))
    else:
        if np.sign(alpha) == 1:
            alpha = alpha % (2 * np.pi)
        else:
            alpha = alpha % (-2 * np.pi)
    #print(np.degrees([alpha, alpha1]))
    q_dynamic = 0.5 * rho * (tail_xdot ** 2 + tail_ydot ** 2) * S_tail

    # Sigma --> Blending Function
    sigma_num = (1 + np.exp(-M*(alpha - alpha0)) + np.exp(M*(alpha + alpha0)))
    sigma_denom = (1 + np.exp(-M*(alpha - alpha0))) * (1 + np.exp(M*(alpha + alpha0)))
    sigma = sigma_num / sigma_denom

    CL_0 = 0.0 #0.23
    CD_p = 0.0
    CL_alpha = np.pi*AR_tail/(1+np.sqrt(1 + (AR_tail/2)**2))
    CD_alpha = ((CL_0 + CL_alpha*alpha)**2)/(np.pi*e*AR_tail)


    CL = ((1 - sigma)*(CL_0 + CL_alpha*alpha)) + ((sigma)*2*np.sign(alpha)*(np.sin(alpha)**2)*np.cos(alpha))
    CD = CD_p + CD_alpha
    #CD = 1-np.cos(2*alpha)

    #CL = np.sin(2*alpha) #2 * np.pi * alpha
    #CD = 1-np.cos(2*alpha) #np.abs(1.28 * np.sin(alpha))
    f_lift = CL * q_dynamic
    f_drag = CD * q_dynamic
    fx_tail = (-f_drag * np.cos(alpha) + (f_lift * np.sin(alpha)))
    fy_tail = (f_drag * np.sin(alpha) + (f_lift * np.cos(alpha)))
    # fx_tail = (-f_drag * np.cos(alpha) + -(f_lift * np.sin(alpha)))
    # fy_tail = (-f_drag * np.sin(alpha) + (f_lift * np.cos(alpha)))
    fx = fx_tail*np.cos(angle_tail) + fy_tail*np.sin(angle_tail)
    fy = -fx_tail*np.sin(angle_tail) + fy_tail*np.cos(angle_tail)
    #print(f"Alpha: {np.degrees(alpha)}, Lift: {f_lift}, Drag: {f_drag}, fx_tail: {fx_tail}, fy_tail: {fy_tail}, fx: {fx}, fy: {fy}")


    # Handle the case when inputs were scalars
    if isinstance(t, (int, float)):
        return fx.item(), fy.item(), alpha.item()  # Convert scalars back to Python scalars
    else:
        return np.array([fx, fy, alpha])  # Return NumPy arrays for arrays


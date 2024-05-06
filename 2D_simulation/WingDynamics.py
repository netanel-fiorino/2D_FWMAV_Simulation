import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('..')

from parameters import *
from WingShape import WingShape

class WingDynamics:
    def __init__(self, side, Wingshape, num_elements=num_elements, glide=False, time_step=t_step, \
                 v0=v_free_stream, freq=freq, \
                 phi_amplitude=phi_amplitude, theta_amplitude=theta_amplitude, \
                 beta=beta, phi=phi_start, psi=psi_start, theta=theta_start, \
                 phi_dot=0, psi_dot=0, theta_dot=0, \
                 phi_0=phi0, psi0=psi0, theta0=theta0):
        self._state = np.array([[0],    # (x-pos) [0]
                               [0],     # (y-pos) [1]
                               [0],     # (z-pos) [2]
                               [v0.item(0)],     # (u) [3]
                               [v0.item(1)],     # (v) [4]
                               [v0.item(2)],     # (w) [5]
                               [beta],     # (beta) angle to stroke plane [6]
                               [phi],     # (phi) causes flapping motion [7]
                               [psi],     # (psi) deviation angle [8]
                               [theta],      # (theta) rotation angle[9]
                               [phi_dot],      # (phi_dot) [10]
                               [psi_dot],      # (psi_dot) [11]
                               [theta_dot]])     # (theta_dot) [12]
        self.side = side
        self.Wingshape = Wingshape
        self.time_step = time_step
        self.num_elements = num_elements

        self.freq = freq
        self.phi_amplitude = phi_amplitude
        self.theta_amplitude = theta_amplitude
        self.phi0 = phi0         #mean stroke angle
        self.theta0 = theta0     #mean rotation angle
        self.psi0 = psi0     #mean deviation angle


        self._state[10][0] = -self.phi_amplitude*2*np.pi*self.freq*np.cos(2*np.pi*self.freq*0 + np.pi/2 )                                                                       #derivative of phi function at t=t_period
        self._state[12][0] = 2*np.pi*self.theta_amplitude*C_alpha*self.freq * (np.tanh(C_alpha))**(-1) * np.cos(2*np.pi*self.freq*0) * (np.cosh(-C_alpha*np.sin(2*np.pi*self.freq*0))**-2)     #derivative of theta function at t=t_period

        self.u_inf = np.sqrt(self._state.item(3)**2 + self._state.item(4)**2 + self._state.item(5)**2)
        self.J = self.u_inf/(2*(2*self.phi_amplitude)*self.Wingshape.R*self.freq)

        self.V_W_inflow_prev = np.zeros((self.num_elements, 3, 1))

        self.glide = glide
        if self.glide == True:
            self._state[10][0] = 0
            self._state[11][0] = 0
            self._state[12][0] = 0


    def update(self, time):
        if self.glide == False:
            w = 2*np.pi*self.freq

            new_stroke = self.phi0 - self.phi_amplitude*np.sin(2*np.pi*self.freq*time + np.pi/2)#-np.pi/2)
            #new_stroke = 0.1192  + -0.3018*np.cos(w*time) + -0.3905*np.sin(w*time) + -0.01042*np.cos(2*w*time) + -0.02967*np.sin(2*w*time) + -0.002753*np.cos(3*w*time) + -0.001312*np.sin(3*w*time)
            #https://www.mdpi.com/2226-4310/3/3/23

            C_alpha = 2.6
            new_pitch = self.theta0 - self.theta_amplitude/(np.tanh(C_alpha))*np.tanh(C_alpha*np.sin(2*np.pi*self.freq*time -np.pi/2 +np.pi -np.pi/2 + np.pi))           #taken from Bhatti
            #new_pitch = -0.09502  + -0.5728*np.cos(w*time) + -0.479*np.sin(w*time) + -0.01803*np.cos(2*w*time) + -0.0364*np.sin(2*w*time) + -0.1374*np.cos(3*w*time) + -0.1374*np.sin(3*w*time)

            new_deviation = self.psi0
            self._state[10][0] = (new_stroke - self._state[7][0])/self.time_step
            self._state[11][0] = (new_deviation - self._state[8][0])/self.time_step
            self._state[12][0] = (new_pitch - self._state[9][0])/self.time_step

            self._state[7][0] = new_stroke
            self._state[8][0] = new_deviation
            self._state[9][0] = new_pitch
        else:
            self._state[10][0] = 0
            self._state[11][0] = 0
            self._state[12][0] = 0

            self._state[7][0] = self._state.item(7)
            self._state[8][0] = self._state.item(8)
            self._state[9][0] = self._state.item(9)

    def update_initial_angle(self, time, stroke_offset, pitch_offset):
        if self.glide == False:
            w = 2*np.pi*self.freq

            new_stroke = self.phi0 - self.phi_amplitude*np.sin(2*np.pi*self.freq*time + np.pi/2 - stroke_offset)
            #new_stroke = 0.1192  + -0.3018*np.cos(w*time) + -0.3905*np.sin(w*time) + -0.01042*np.cos(2*w*time) + -0.02967*np.sin(2*w*time) + -0.002753*np.cos(3*w*time) + -0.001312*np.sin(3*w*time)
            #https://www.mdpi.com/2226-4310/3/3/23

            C_alpha = 2.6
            new_pitch = self.theta0 - self.theta_amplitude/(np.tanh(C_alpha))*np.tanh(C_alpha*np.sin(2*np.pi*self.freq*time + np.pi - pitch_offset))           #taken from Bhatti
            #new_pitch = -0.09502  + -0.5728*np.cos(w*time) + -0.479*np.sin(w*time) + -0.01803*np.cos(2*w*time) + -0.0364*np.sin(2*w*time) + -0.1374*np.cos(3*w*time) + -0.1374*np.sin(3*w*time)

            new_deviation = self.psi0
            self._state[10][0] = (new_stroke - self._state[7][0])/self.time_step
            self._state[11][0] = (new_deviation - self._state[8][0])/self.time_step
            self._state[12][0] = (new_pitch - self._state[9][0])/self.time_step

            self._state[7][0] = new_stroke
            self._state[8][0] = new_deviation
            self._state[9][0] = new_pitch
        else:
            self._state[10][0] = 0
            self._state[11][0] = 0
            self._state[12][0] = 0

            self._state[7][0] = self._state.item(7)
            self._state[8][0] = self._state.item(8)
            self._state[9][0] = self._state.item(9)

    def get_aero_coeff(self, r):
        self.u_inf = self._state.item(3)#np.sqrt(self._state.item(3)**2 + self._state.item(4)**2 + self._state.item(5)**2)
        self.J = self.u_inf/(2 * (2*self.phi_amplitude) * self.Wingshape.R * self.freq)
        xhat_0 = self.Wingshape.Y_LE_hat__r_hat(r/self.Wingshape.R) + self.Wingshape.y_hat__r      #=x/c where x is the distance between the leading edge and the wing pitching axis
        K_PL = -2.109*((self.J + self.Wingshape.rhat2)**-0.606) + 4.136
        K_VL = 2.659*((self.J + self.Wingshape.rhat2)**-0.666) + -0.344
        K_PD = -0.182*((self.J + self.Wingshape.rhat2)**-2.414) + 1.370
        K_VD = 0.765*((self.J + self.Wingshape.rhat2)**-1.497) + 2.078
        K_PM = 0.803*((self.J + self.Wingshape.rhatM)**-0.972) + -0.363
        K_VM = -0.242*((self.J + self.Wingshape.rhatM)**-1.354) + -0.554

        alpha = self.alpha
        # if alpha>np.radians(90):
        #     alpha = -np.pi + alpha
        
        CL = K_PL*np.sin(alpha)*(np.cos(alpha)**2) + K_VL*(np.sin(alpha)**2)*np.cos(alpha)
        CD = K_PD*(np.sin(alpha)**2)*np.cos(alpha) + K_VD*(np.sin(alpha)**3)
        CM = K_PM*(np.sin(alpha)**2)*np.cos(alpha) + K_VM*(np.sin(alpha)**2)
        CR = np.pi*(0.75-xhat_0)
        CA = np.pi/8
        return CL, CD, CM, CR, CA

    def force_calc(self):
        beta = self._state.item(6)
        phi = self._state.item(7)
        psi = self._state.item(8)
        theta = self._state.item(9)
        phi_dot = self._state.item(10)
        psi_dot = self._state.item(11)
        theta_dot = self._state.item(12)

        if self.side == 'right':              #side when looking at fwmav, (side along negative y axis)
            #Dont know if all these should be negatives
            phi = -phi
            phi_dot = -phi_dot
            theta = theta
            theta_dot = theta_dot
            psi = -psi
            psi_dot = -psi_dot

        R_BW = R_theta(theta) @ R_psi(psi) @ R_phi(phi) @ R_beta(beta)

        #convention for naming: V_W_body is velocity of the body expressed in the wing frame (W)
        V_W_body = R_BW @ self._state[3:6]      #rotation @ body velocity
        omega_W_wing = np.array([[0], [theta_dot], [0]]) + R_theta(theta) @ np.array([[psi_dot], [0], [0]]) + R_theta(theta) @ R_psi(psi) @ np.array([[0], [0], [phi_dot]])
        if self.glide == False:
            #omega_W_wing = np.array([[0], [theta_dot], [0]]) + R_theta(theta) @ np.array([[psi_dot], [0], [0]]) + R_theta(theta) @ R_psi(psi) @ np.array([[0], [0], [phi_dot]])
            F_W_trans = np.array([[0],[0],[0]])
            F_W_rot = np.array([[0],[0],[0]])
            F_W_added = np.array([[0],[0],[0]])

            M_W_trans = np.array([[0],[0],[0]])
            M_W_rot = np.array([[0],[0],[0]])
            M_W_added = np.array([[0],[0],[0]])

            for i in range(self.num_elements):
                #aero calculations
                delta_r = self.Wingshape.R/self.num_elements
                if self.side == 'left':                                            #side when looking at fwmav, (side along positive y axis)
                    r_W_i = np.array([[0], [(i+0.5)*delta_r], [0]])                 #distance from wing pivot to i-th blade element
                    chord_i = self.Wingshape.C_hat__r_hat(r_W_i.item(1)/self.Wingshape.R)*self.Wingshape.C_mean

                    V_W_inflow = V_W_body + np.cross(omega_W_wing.T, r_W_i.T).T
                    V_W_i = np.array([[1, 0, 0],[0, 0, 0],[0, 0, 1]]) @ V_W_inflow
                    chat_W_i = np.array([[0], [0], [-1]])                   #c hat in Wing frame for element i

                    
                    V_W_COP = np.array([[1, 0, 0],[0, 0, 0],[0, 0, 1]]) @ (V_W_body + np.cross(omega_W_wing.T, np.array([[0], [self.Wingshape.rhat2*self.Wingshape.R], [0]]).T).T)
                    self.alpha_COP = np.arctan2((np.sqrt(np.cross(V_W_COP.T, chat_W_i.T) @ np.cross(V_W_COP.T, chat_W_i.T).T)),(V_W_COP.T @ chat_W_i))

                    self.alpha = np.arctan2((np.sqrt(np.cross(V_W_i.T, chat_W_i.T) @ np.cross(V_W_i.T, chat_W_i.T).T)),(V_W_i.T @ chat_W_i))
                    
                    CL, CD, CM, CR, CA = self.get_aero_coeff(r_W_i.item(1))
                    #calculate frames
                    i_W = np.array([[1], [0], [0]])
                    j_W = np.array([[0], [1], [0]])
                    c_i = (self.Wingshape.Y_LE_hat__r_hat(r_W_i.item(1)/self.Wingshape.R) - self.Wingshape.C_hat__r_hat(r_W_i.item(1)/self.Wingshape.R)/2) * self.Wingshape.C_mean          #distance from wing root axis to COM of element
                    s_W_i = np.array([[0], [np.dot(r_W_i.T, j_W).item(0)], [-c_i]])                                                                                                          #Location of center of blade element
                    lhat_W_i = (np.dot(V_W_i.T, i_W))/(np.sqrt(np.dot(V_W_i.T, i_W) * np.dot(V_W_i.T, i_W))) * np.array([[0, 0, 1], [0, 1, 0], [-1, 0, 0]]) @ (V_W_i)/(np.sqrt(V_W_i.T @ V_W_i))
                    dhat_W_i = np.array([[-1, 0, 0], [0, 0, 0], [0, 0, -1]]) @ (V_W_i)/(np.sqrt(V_W_i.T @ V_W_i))

                    V_W_inflow_s = V_W_body + np.cross(omega_W_wing.T, s_W_i.T).T
                    a_W_i_inflow = (V_W_inflow_s - self.V_W_inflow_prev[i])/self.time_step        #acceleration of center of blade element (derivative of velocity)
                    self.V_W_inflow_prev[i,:] = V_W_inflow_s

                    F_W_trans_i = (CL* 0.5 * rho * (V_W_i.T @ V_W_i) * chord_i * delta_r) * lhat_W_i + (CD * 0.5 * rho * (V_W_i.T @ V_W_i) * chord_i * delta_r) * dhat_W_i
                    F_W_rot_i = (CR * rho * np.dot(omega_W_wing.T, j_W) * np.sqrt(V_W_i.T @ V_W_i) * chord_i**2 * delta_r) * i_W
                    F_W_added_i = (CA * rho * np.dot(a_W_i_inflow.T, i_W) * chord_i**2 * delta_r) * i_W

                    F_W_trans = F_W_trans + F_W_trans_i
                    F_W_rot = F_W_rot + F_W_rot_i
                    F_W_added = F_W_added + F_W_added_i
                    
                    M_W_trans_i = (CM * 0.5 * rho * V_W_i.T @ V_W_i * chord_i**2 * delta_r) * j_W + np.cross(r_W_i.T, F_W_trans_i.T).T
                    M_W_rot_i = np.cross(s_W_i.T, F_W_rot_i.T).T
                    M_W_added_i = np.cross(s_W_i.T, F_W_added_i.T).T

                    
                    M_W_trans = M_W_trans + M_W_trans_i
                    M_W_rot = M_W_rot + M_W_rot_i
                    M_W_added = M_W_added + M_W_added_i

                elif self.side == 'right':
                    r_W_i = np.array([[0], [(i+0.5)*delta_r], [0]])
                    chord_i = self.Wingshape.C_hat__r_hat(r_W_i.item(1)/self.Wingshape.R)*self.Wingshape.C_mean

                    V_W_inflow = V_W_body + np.cross(omega_W_wing.T, -r_W_i.T).T
                    V_W_i = np.array([[1, 0, 0],[0, 0, 0],[0, 0, 1]]) @ V_W_inflow
                    chat_W_i = np.array([[0], [0], [-1]])                   #c hat in Wing frame for element i

                    
                    V_W_COP = np.array([[1, 0, 0],[0, 0, 0],[0, 0, 1]]) @ (V_W_body + np.cross(omega_W_wing.T, -np.array([[0], [self.Wingshape.rhat2*self.Wingshape.R], [0]]).T).T)
                    self.alpha_COP = np.arctan2((np.sqrt(np.cross(V_W_COP.T, chat_W_i.T) @ np.cross(V_W_COP.T, chat_W_i.T).T)),(V_W_COP.T @ chat_W_i))

                    self.alpha = np.arctan2((np.sqrt(np.cross(V_W_i.T, chat_W_i.T) @ np.cross(V_W_i.T, chat_W_i.T).T)),(V_W_i.T @ chat_W_i))
                    #self.alpha = np.arctan((np.sqrt(np.cross(V_W_i.T, chat_W_i.T) @ np.cross(V_W_i.T, chat_W_i.T).T))/(V_W_i.T @ chat_W_i))
                    
                    CL, CD, CM, CR, CA = self.get_aero_coeff(r_W_i.item(1))
                    
                    #calculate frames
                    i_W = np.array([[1], [0], [0]])
                    j_W = np.array([[0], [1], [0]])
                    c_i = (self.Wingshape.Y_LE_hat__r_hat(r_W_i.item(1)/self.Wingshape.R) - self.Wingshape.C_hat__r_hat(r_W_i.item(1)/self.Wingshape.R)/2) * self.Wingshape.C_mean          #distance from pitching axis to COM of element
                    s_W_i = np.array([[0], [np.dot(-r_W_i.T, j_W).item(0)], [-c_i]])
                    lhat_W_i = (np.dot(V_W_i.T, i_W))/(np.sqrt(np.dot(V_W_i.T, i_W) * np.dot(V_W_i.T, i_W))) * np.array([[0, 0, 1], [0, -1, 0], [-1, 0, 0]]) @ (V_W_i)/(np.sqrt(V_W_i.T @ V_W_i))
                    dhat_W_i = np.array([[-1, 0, 0], [0, 0, 0], [0, 0, -1]]) @ (V_W_i)/(np.sqrt(V_W_i.T @ V_W_i))

                    V_W_inflow_s = V_W_body + np.cross(omega_W_wing.T, s_W_i.T).T
                    a_W_i_inflow = (V_W_inflow_s - self.V_W_inflow_prev[i])/self.time_step        #acceleration of center of blade element (derivative of velocity)
                    self.V_W_inflow_prev[i,:] = V_W_inflow_s

                    F_W_trans_i = (CL* 0.5 * rho * (V_W_i.T @ V_W_i) * chord_i * delta_r) * lhat_W_i + (CD * 0.5 * rho * (V_W_i.T @ V_W_i) * chord_i * delta_r) * dhat_W_i
                    F_W_rot_i = (CR * rho * np.dot(omega_W_wing.T, j_W) * np.sqrt(V_W_i.T @ V_W_i) * chord_i**2 * delta_r) * i_W
                    F_W_added_i = (CA * rho * np.dot(a_W_i_inflow.T, i_W) * chord_i**2 * delta_r) * i_W

                    F_W_trans = F_W_trans + F_W_trans_i
                    F_W_rot = F_W_rot + F_W_rot_i
                    F_W_added = F_W_added + F_W_added_i

                    
                    M_W_trans_i = (CM * 0.5 * rho * V_W_i.T @ V_W_i * chord_i**2 * delta_r) * j_W + np.cross(-r_W_i.T, F_W_trans_i.T).T
                    M_W_rot_i = np.cross(s_W_i.T, F_W_rot_i.T).T
                    M_W_added_i = np.cross(s_W_i.T, F_W_added_i.T).T

                    M_W_trans = M_W_trans + M_W_trans_i
                    M_W_rot = M_W_rot + M_W_rot_i
                    M_W_added = M_W_added + M_W_added_i

            #print(F_W_added, F_W_rot, F_W_trans)
            F_W = F_W_trans + F_W_rot + F_W_added
            M_W = M_W_trans + M_W_rot + M_W_added

            F_B = R_BW.T @ F_W
            M_B = R_BW.T @ M_W


            return np.array([[F_B.item(0), F_B.item(1), F_B.item(2), M_B.item(0), M_B.item(1), M_B.item(2)]]).T
        elif self.glide == True:
            chat_W_i = np.array([[0], [0], [-1]])                   #c hat in Wing frame for element i
            #omega_W_wing = np.array([[0], [theta_dot], [0]]) + R_theta(theta) @ np.array([[psi_dot], [0], [0]]) + R_theta(theta) @ R_psi(psi) @ np.array([[0], [0], [phi_dot]])

            delta_r = self.Wingshape.R
            r_W_i = np.array([[0], [self.Wingshape.rhat2*delta_r], [0]])                 #distance from wing pivot to i-th blade element
            #chord_i = self.Wingshape.C_hat__r_hat(r_W_i.item(1))*self.Wingshape.C_mean
            if self.side == 'left':                                            #side when looking at fwmav, (side along positive y axis)
                 V_W_inflow = V_W_body + np.cross(omega_W_wing.T, r_W_i.T).T                 #distance from wing pivot to i-th blade element
            elif self.side == 'right':
                 V_W_inflow = V_W_body + np.cross(omega_W_wing.T, -r_W_i.T).T
            V_W_i = np.array([[1, 0, 0],[0, 0, 0],[0, 0, 1]]) @ V_W_inflow
            #print(V_W_i)
            if (V_W_i.T @ chat_W_i) == 0:
                self.alpha_COP = np.pi/2
            else:
                self.alpha_COP = np.arctan((np.sqrt(np.cross(V_W_i.T, chat_W_i.T) @ np.cross(V_W_i.T, chat_W_i.T).T))/(V_W_i.T @ chat_W_i))
                #self.alpha_COP = np.arctan2((np.sqrt(np.cross(V_W_i.T, chat_W_i.T) @ np.cross(V_W_i.T, chat_W_i.T).T)),(V_W_i.T @ chat_W_i))
            self.alpha_COP = np.arctan2(V_W_i.item(0),-V_W_i.item(2))

            q_dynamic = 0.5 * rho * (V_W_body[0] ** 2 + V_W_body[2] ** 2) * self.Wingshape.S

            # Sigma --> Blending Function
            sigma_num = (1 + np.exp(-M*(self.alpha_COP - alpha0)) + np.exp(M*(self.alpha_COP + alpha0)))
            sigma_denom = (1 + np.exp(-M*(self.alpha_COP - alpha0))) * (1 + np.exp(M*(self.alpha_COP + alpha0)))
            sigma = sigma_num / sigma_denom

            CL_0 = 0 #0.23
            CL_alpha = np.pi*self.Wingshape.AR/(1+np.sqrt(1 + (self.Wingshape.AR/2)**2))
            CD_alpha = ((CL_0 + CL_alpha*self.alpha_COP)**2)/(np.pi*e*self.Wingshape.AR)

            CL = ((1 - sigma)*(CL_0 + CL_alpha*self.alpha_COP)) + ((sigma)*2*np.sign(self.alpha_COP)*(np.sin(self.alpha_COP)**2)*np.cos(self.alpha_COP))
            CD = CD_alpha

            f_lift = CL * q_dynamic
            f_drag = CD * q_dynamic
            fz = f_drag * np.cos(self.alpha_COP) + (-f_lift * np.sin(self.alpha_COP))
            fx = (-f_drag * np.sin(self.alpha_COP) + (-f_lift * np.cos(self.alpha_COP)))

            F_W = np.array([[fx.item(0)],[0],[fz.item(0)]])
            M_W = np.array([[0],[0],[0]])
            F_B = R_BW.T @ F_W
            M_B = R_BW.T @ M_W


            return np.array([[F_B.item(0), F_B.item(1), F_B.item(2), M_B.item(0), M_B.item(1), M_B.item(2)]]).T



def R_theta(theta):
    c_theta = np.cos(theta)
    s_theta = np.sin(theta)
    return np.array([[c_theta, 0, -s_theta],
                        [0, 1, 0],
                        [s_theta, 0, c_theta]])
def R_psi(psi):
    c_psi = np.cos(psi)
    s_psi = np.sin(psi)
    return np.array([[1, 0, 0],
                      [0, c_psi, -s_psi],
                      [0, s_psi, c_psi]])
def R_phi(phi):
    c_phi = np.cos(phi)
    s_phi = np.sin(phi)
    return np.array([[c_phi, s_phi, 0],
                      [-s_phi, c_phi, 0],
                      [0, 0, 1]])
def R_beta(beta):
    c_beta = np.cos(beta)
    s_beta = np.sin(beta)
    return np.array([[c_beta, 0, s_beta],
                        [0, 1, 0],
                        [-s_beta, 0, c_beta]])
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('..')

from parameters import *

class WingShape:
    def __init__(self, R, AR, y_hat__r):
        #R is the projected distance along the r'-axis from the wing root to the most distal point on the wing
        #AR is the aspect ratio of the wing
        #r_hat is the nondimensional radial distance defined as r'/R
        #C_hat is the nondimensional chord profile defined as chord/chord_mean
        #C_hat__r_hat is the non dimensional chord profile as a function of the nondimensional radial distance r_hat
        #Y_LE_hat is the non dimensional leading edge profile defined as Y_LE/chord_mean
        #Y_LE_hat__r_hat is the nondimensional leading edge profile as a function of the nondimensional radial distance r_hat

        self.R = R                          #projected distance along the r'-axis from the wing root to the most distal point on the wing
        self.AR = AR                       #Aspect Ratio of the wing
        self.C_mean = self.R / self.AR  #area of one wing divided by wing length
        self.rhat2()
        self.rhatM()
        self.S = self.R * self.C_mean
        self.y_hat__r = y_hat__r



    def chord_f(self, r):
        if self.wing_shape == 'rectangle':
            return self.max_chord                                            #rectangle wing
        elif self.wing_shape == 'ellipse':
            return self.max_chord*np.sqrt(1-r**2/(self.wing_length**2))           #ellipse wing


    def C_hat__r_hat(self, r_hat):
        C_hat = 1               #for rectangular wing
        #C_hat = np.abs((-0.01305/self.C_mean-(-0.09525/self.C_mean))*r_hat - 0.09525/self.C_mean)          #(end_chord_length/C_mean-front_chord_length/C_mean)*r_hat - front_chord_length/C_mean
        #C_hat = np.abs(2*r_hat - 2)
        #C_hat = np.abs(1*r_hat - 1.5)
        return C_hat

    def Y_LE_hat__r_hat(self, r_hat):
        Y_LE_hat = 0.0
        return Y_LE_hat

    def rhat2(self):
        rhat2 = 0
        for i in range(1000):
            deltar = 1/1000
            rhat2 = rhat2 + ((i+0.5)*deltar)**2 * self.C_hat__r_hat((i+0.5)*deltar) * deltar
        self.rhat2 = np.sqrt(rhat2)

    def rhatM(self):
        rhatM = 0
        for i in range(1000):
            deltar = 1/1000
            rhatM = rhatM + ((i+1)*deltar)**2 * (self.C_hat__r_hat((i+1)*deltar))**2 * deltar
        self.rhatM = np.sqrt(rhatM)
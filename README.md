# 2D_FWMAV_Simulation
FWMAV Physics Simulation

## Overview
This repository contains the code developed as part of my Master's thesis project at The Cooper Union for the Advancement of Science and Art. The project focuses on the development of a physics simulation for a Flapping Wing Micro Aerial Vehicle (FWMAV).

## Usage
To use the simulation software, follow these steps:

1. Parameters: The parameters for the simulation can be found in the parameter file. You can modify these parameters to adjust for various FWMAVs and change the simulation time.
2. Wing Shape: Changes to the wing shape can be made in the WingShape class. Modify the parameters in this class to experiment with different geometries.
3. Wing Kinematics: Modify the equations governing wing kinematics in the 'update_initial_angle' function within the WingDynamics class. This function updates the wing kinematics based on the initial wing angles.
4. Initial Flight Conditions: Set the initial flight conditions in the main file. Adjust parameters such as initial position, velocity, and orientation to simulate different flight scenarios.
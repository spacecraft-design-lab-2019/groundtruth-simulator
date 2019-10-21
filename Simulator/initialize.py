# -*- coding: utf-8 -*-
"""
Initialization Script

Sets simulation parameters and initializes spacecraft. Run by driver-script
before beginning simulation.

"""

#-------------------------Setup---------------------------------

# Importing libraries
import numpy as np

# Simulation Parameters
tspan = np.array([0, 864])    # [sec]
tstep = .1                     # [sec] - 10 Hz


# Initial Spacecraft State
r_i = np.array([6712880.93e-3,1038555.54e-3,-132667.04e-3])
q_i = np.array([1, 0, 0, 0])
v_i = np.array([-831.937369e-3,4688.525767e-3,-6004.570270e-3])
w_i = np.array([.1, .5, -.3])
state_i = np.r_[r_i, q_i, v_i, w_i]


#--------------------------------------------------------------



# class SpacecraftState():
#     """
#     A class to store spacecraft position, velocity, attitude, and angular velocity.
#     """
#     def __init__(self, r=r_i, q=q_i, v=v_i, w=w_i):
#         self.state = np.r_[r, q, v, w]

#     def r(self):
#         return state[0:3]
#     def q(self):
#         return state[3:7]
#     def v(self):
#         return state[7:10]
#     def w(self):
#         return state[10:13]

#     def normalize_quat(self):
#         self.state[3:7] = self.q()/np.linalg.norm(self.q())

#     def __add__(self, other):
#         sum = self + other
#         normalize_quat()


# class SpacecraftState():
#     """
#     A class to store all spacecraft parameters
#     """
#     def __init__(self,
#                 I = np.array([[17,0,0],[0,18,0],[0,0,22]]),
#                 r = np.array([6712880.93e-3,1038555.54e-3,-132667.04e-3]),
#                 q = np.array([1, 0, 0, 0]), # scalar first
#                 v = np.array([-831.937369e-3,4688.525767e-3,-6004.570270e-3]),
#                 w = np.array([0, 0, 0]),
#                 t = 0): #time
#         self.I = I
#         self.r = r
#         self.q = q
#         self.v = v
#         self.w = w
#         self.t = t

#     def state(self):
#         return np.r_[self.r, self.q, self.v, self.w]

#     def set_state(self,new_state, t):
#         self.r = new_state[0:3]
#         self.q = new_state[3:7]
#         self.v = new_state[7:10]
#         self.w = new_state[10:13]
#         self.t = t
#         return self

#     def update_state(self, state_dot, dt):
#         self.r = self.r + state_dot[0:3]
#         self.q = self.q + state_dot[3:7]
#         self.v = self.v + state_dot[7:10]
#         self.w = self.w + state_dot[10:13]
#         self.t = self.t + dt
#         return self

#     def propogate_state(self, dself, dt):
#         x = copy(self)
#         x.update_state(self, self+dself, self.tdt)

#     def normalize_quat(self):
#         self.q = self.q/np.linalg.norm(self.q)

#     def __add__(self, other):


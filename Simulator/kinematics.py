# -*- coding: utf-8 -*-
import numpy as np
from conversions import skew
#---------------------------Kinematics---------------------------------

def calc_q_dot(q, w):
    """
    Function: calc_q_dot
        Calculates the derivative of the quaternion vector.

    Inputs:
        q: quaternion vector np.array[4x1][float64]
        w: angular velocity (radians/s) np.array[3x1][float64]
    Outputs:
        q_dot: rate of change of quaternion np.array[4x1][float64]
    """
    omega = np.array([[0, w[2], -w[1], w[0]],
                    [-w[2], 0, w[0], w[1]],
                    [w[1], -w[0], 0, w[2]],
                    [-w[0], -w[1], -w[2],0]])
    return 0.5* omega @ q


def calc_w_dot(w, torque, I):
    """
    Function: calc_w_dot
        Calculates the derivative of the angular velocity vector using the
        Euler equations for rigid body dynamics.

    Inputs:
        w: angular velocity (radians/s) np.array[3x1][float64]
        torque: applied torque (Nm) np.array[3x1][float54]
    Ouputs:
        w_dot: derivative of angular velocity (rad/s/s) np.array[3x1][float64]
    """

    # w_dot = np.linalg.solve(I, -np.dot(skew(w), np.dot(I, w)) + torque)
    w_dot = np.linalg.inv(I) @ -(skew(w) @ (I @ w)) + torque
    return w_dot


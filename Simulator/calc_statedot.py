# -*- coding: utf-8 -*-
import numpy as np

def calc_statedot(t, s):
    """
    Function: state_dot
        Calculates the derivative of the state vector
        
    Inputs:
        t: the current time
        s: the current state vector
    Outputs:
        state_dot: the derivative of the state vector
    """
    
    # this is a dummy 
    con = -mu/(np.sqrt(s[0]**2 + s[1]**2 + s[2]**2)**3)
    return np.array([s[3],s[4],s[5],con*s[0],con*s[1],con*s[2]])

    # TO IMPLEMENT:
    # Calculate environment given current JD and state
    
    # Calculate Accelerations/Torques
    
    # Build statedot
    statedot = np.zeros(np.shape(s))
    statedot[0:3] = s[7:10]
    return statedot
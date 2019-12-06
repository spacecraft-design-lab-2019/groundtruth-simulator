import numpy as np
from numpy import linalg as LA
import math

def SkewSymmetric(a):
    a1 = a[0]
    a2 = a[1]
    a3 = a[2]
    a_skew = np.array([[0,-a3,a2],[a3,0,-a1],[-a2,a1,0]])
    return a_skew

def predict(xk,w,dt):
    '''
    Predict step of MEKF - state x = [q, Beta]
    qk*[cos(theta),rsin(theta)] uses Hamilton's quaternion mult.
    ** scalar first quaternion ** 
    '''
    q = np.reshape(xk[0:4],(4,1))
    b = np.reshape(xk[4:7],(3,1))
    w = np.reshape(w,(3,1))
    u = w + b
    theta = LA.norm(u-b)*dt
    r = (u-b)/LA.norm(u-b)
    xn =np.zeros(7)
    
    # q2 == s
    s_scalar = math.cos(theta/2)
    s_vector = r*math.sin(theta/2)
    s_vector = np.reshape(s_vector,(3,1))
    s = np.vstack((s_scalar,s_vector))

    # propagate dynamics forward
    xn[0] = q[0]*s[0]-np.matmul(q[1:4].T,s[1:4])
    temp = q[0]*s[1:4]+s[0]*q[1:4]+np.reshape(np.cross(q[1:4].T,s[1:4].T),(3,1))
    xn[1:4] = np.reshape(temp,3)
    xn[4:7] = np.reshape(b,3)
    # isolate vector component of quaternion
    V = np.array([[0,1,0,0],[0,0,1,0],[0,0,0,1]])
    
    # left quaternion multiply
    Ls = np.zeros([4,4])
    Ls[0,0] = s[0]
    Ls[0,1:4] = -np.reshape(s[1:4],3)
    Ls[1:4,0] = np.reshape(s[1:4],3)
    Ls[1:4,1:4] = s[0]*np.identity(3) + SkewSymmetric(s[1:4])
    
    # right quaternion multiply
    Rs = np.zeros([4,4])
    Rs[0,0] = s[0]
    Rs[0,1:4] = -np.reshape(s[1:4],3)
    Rs[1:4,0] = np.reshape(s[1:4],3)
    Rs[1:4,1:4] = s[0]*np.identity(3) - SkewSymmetric(s[1:4])
    
    
    A = np.zeros([6,6])
    # fill out A matrix
    
    A[0:3,0:3] = V@Ls.T@Rs@V.T
    A[0:3,3:6] = -0.5*dt*np.identity(3)
    A[3:6,3:6] = np.identity(3)

    
        
    return xn, A




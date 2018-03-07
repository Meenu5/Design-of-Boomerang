import numpy as np
from numpy.linalg import inv

tan = np.tan
cos = np.cos
sin = np.sin
pi = np.pi
transpose = np.transpose

# From initial frame to body frame
def doT0Transformation(phi, theta, psi) :
    T0 = np.zeros([3,3])
    T0[0,0] = cos(theta)*cos(psi)
    T0[0,1] = cos(theta)*sin(psi)
    T0[0,2] = -sin(theta)
    T0[1,0] = -cos(phi)*sin(psi)+sin(phi)*sin(theta)*cos(psi)
    T0[1,1] = cos(phi)*cos(psi)+sin(phi)*sin(theta)*sin(psi)
    T0[1,2] = sin(phi)*cos(theta)
    T0[2,0] = sin(phi)*sin(psi)+cos(phi)*sin(theta)*cos(psi)
    T0[2,1] = -sin(phi)*cos(psi)+cos(phi)*sin(theta)*sin(psi)
    T0[2,2] = cos(phi)*cos(theta)
    return T0

# From body frame to non spinning frame
def doTnTransformation(lamda) :
    Tn = np.zeros([3,3])
    Tn[0,0] = cos(lamda)
    Tn[0,1] = -sin(lamda)
    Tn[0,2] = 0
    Tn[1,0] = sin(lamda)
    Tn[1,1] = cos(lamda)
    Tn[1,2] = 0
    Tn[2,0] = 0
    Tn[2,1] = 0
    Tn[2,2] = 1
    return Tn

# From initial frame to non spinning frame
def doT0TnTransformation(T0,Tn) :
    return np.matmul(T0,Tn)

# Transformation from body to blade frame
def doT0TjTransformation(Lamda,theta_pitch,beta) :
    Tj = np.zeros([3,3])
    Tj[0,0] = cos(theta_pitch)*sin(Lamda)+sin(theta_pitch)*sin(beta)*cos(Lamda)
    Tj[0,1] = -cos(theta_pitch)*cos(Lamda)+sin(theta_pitch)*sin(beta)*sin(Lamda)
    Tj[0,2] = -sin(theta_pitch)*cos(beta)
    Tj[1,0] = -cos(beta)*cos(Lamda)
    Tj[1,1] = cos(beta)*sin(Lamda)
    Tj[1,2] = sin(beta)
    Tj[2,0] = sin(theta_pitch)*sin(Lamda)-cos(theta_pitch)*sin(beta)*cos(Lamda)
    Tj[2,1] = -sin(theta_pitch)*cos(Lamda)-cos(theta_pitch)*sin(beta)*sin(Lamda)
    Tj[2,2] = cos(theta_pitch)*cos(beta)
    return Tj

# Transformation from intertial frame to non spinning frame
def doTiTransformation(Phi, Theta, Psi) :
    Ti = np.zeros([3,3])
    Ti[0,0] = cos(Theta)*cos(Psi)
    Ti[0,1] = cos(Theta)*sin(Psi)
    Ti[0,2] = -sin(Theta)
    Ti[1,0] = -cos(Phi)*sin(Psi)+sin(Phi)*sin(Theta)*cos(Psi)
    Ti[1,1] = cos(Phi)*cos(Psi)+sin(Phi)*sin(Theta)*sin(Psi)
    Ti[1,2] = sin(Phi)*cos(Theta)
    Ti[2,0] = sin(Phi)*sin(Psi)+cos(Phi)*sin(Theta)*cos(Psi)
    Ti[2,1] = -sin(Phi)*cos(Psi)+cos(Phi)*sin(Theta)*sin(Psi)
    Ti[2,2] = cos(Phi)*cos(Theta)
    return Ti

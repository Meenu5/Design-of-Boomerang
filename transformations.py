import numpy as np
from numpy.linalg import inv

tan = np.tan
cos = np.cos
sin = np.sin
pi = np.pi
transpose = np.transpose
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

def doT0TnTransformation(T0,Tn) :
    return np.matmul(T0,Tn)

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

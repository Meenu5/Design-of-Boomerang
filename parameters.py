import numpy as np
import transformations
from numpy.linalg import inv

tan = np.tan
cos = np.cos
sin = np.sin
pi = np.pi
transpose = np.transpose

# f1 - returns segments*3 - position vectors of ac of sections
def doPositionVector(Tj,length,x_ac) :
    r_j_vec = []
    segments = 100
    for i in range(segments) :
        eta = (i+0.5)*length/segments
        pos_vec = np.array([0,eta,0])
        r_j_vec.append(transpose(np.array([x_ac,0,0])) + np.matmul(inv(Tj),transpose(pos_vec)))
    return np.array(r_j_vec)


theta_1 = 0
theta_2 = 0

Lamda_1 = pi*2/3
Lamda_2 = pi*4/3

x_ac = 0.0723

theta_pitch_vec = np.array([theta_1, theta_2])
Lamda_vec = np.array([Lamda_1,Lamda_2])
beta_vec = np.array([0,0])

Tj1 = transformations.doT0TjTransformation(Lamda_vec[0],theta_pitch_vec[0],beta_vec[0]) # blade 1
Tj2 = transformations.doT0TjTransformation(Lamda_vec[1],theta_pitch_vec[1],beta_vec[1]) # blade 2

r_j_vec_1 = []
r_j_vec_2 = []

import numpy as np
from numpy.linalg import inv

tan = np.tan
cos = np.cos
sin = np.sin
pi = np.pi
taninv = np.arctan # first arg/ second arg 's inverse'
transpose = np.transpose

# change in euler angles of body frame
def doEulerRatesBody(p,q,r,phi, theta, psi) :
    psi_d = q*sin(phi)/cos(theta)+r*cos(phi)/cos(theta)
    theta_d = q*cos(phi) - r*sin(phi)
    phi_d = p+ q*sin(phi)*tan(theta)+r*cos(phi)*tan(theta)
    return phi_d, theta_d, psi_d



# f1 - returns segments*3 - position vectors of ac of sections
def doPositionVector(Tj,length,x_ac,segments) :
    r_j_vec = []
    for i in range(segments) :
        eta = (i+0.5)*length/segments
        pos_vec = np.array([0,eta,0])
        r_j_vec.append(transpose(np.array([x_ac,0,0])) + np.matmul(inv(Tj),transpose(pos_vec)))
    return np.array(r_j_vec)

# f6 - Relative air velocity of blade in blade frame
def doRelativeAirVelBlade(v_j_vec, Tj) :
    return transpose(np.matmul(Tj,transpose(-v_j_vec-[0,0,-4.648])))

# f5 - Velocity of blade in body frame
def doVelBlade(u_vec, omega_vec, r_j_ac) :
    v_vec = []
    for i in r_j_ac :
        v_vec.append(u_vec + np.cross(omega_vec,i))
    return np.array(v_vec)

# f7 - Calculation of alpha
def doAlpha(w_vec) :
    alpha_vec = []
    k_vec = []
    for i in w_vec :
        if i[0] < 0 :
            k_vec += [1]
        else :
            k_vec += [-1]

    for i in w_vec :
        if i[2]>=0 and i[0]>=0 :
            alpha = pi - taninv(abs(i[2]/i[0]))
        elif i[2]<=0 and i[0]>=0 :
            alpha = taninv(abs(i[2]/i[0])) - pi
        elif i[2]<=0 and i[0]<=0 :
            alpha = -taninv(abs(i[2]/i[0]))
        else :
            alpha = taninv(abs(i[2]/i[0]))
        alpha_vec.append(alpha)
    return list(alpha_vec),k_vec


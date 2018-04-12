import numpy as np
from numpy.linalg import inv
import aerocoefficients

tan = np.tan
cos = np.cos
sin = np.sin
pi = np.pi
transpose = np.transpose

def doClCd(alpha) :
    Cl, alpha_Cl, Cd, alpha_Cd = aerocoefficients.doLoadLDCoefficients()
    num_Cl = Cl.shape
    num_Cd = Cd.shape

    for i in range(num_Cl[0]) :

        if (abs(alpha_Cl[i]-alpha) < 10**-2) :
            Cl_ans = Cl[i]
            break

        if alpha_Cl[i] > alpha :
            Cl_ans = (Cl[i]-Cl[i-1]) / (alpha_Cl[i]-alpha_Cl[i-1]) * (alpha-alpha_Cl[i-1]) + Cl[i-1]
            break

    for i in range(num_Cd[0]) :

        if (abs(alpha_Cd[i]-alpha) < 10**-4) :
            Cd_ans = Cd[i]
            break

        if alpha_Cd[i] > alpha :
            Cd_ans = (Cd[i]-Cd[i-1]) / (alpha_Cd[i]-alpha_Cd[i-1]) * (alpha-alpha_Cd[i-1]) + Cd[i-1]
            break

    # return 0,0
    return Cl_ans, Cd_ans

def doCm(alpha) :
    Cm, alpha_Cm = aerocoefficients.doLoadMomentCoefficients()
    num_Cm = Cm.shape

    for i in range(num_Cm[0]) :

        if abs(alpha_Cm[i]-alpha) < 10**-4 :
            Cm_ans = Cm[i]
            break

        if alpha_Cm[i] > alpha :
            Cm_ans = (Cm[i]-Cm[i-1]) / (alpha_Cm[i]-alpha_Cm[i-1]) * (alpha-alpha_Cm[i-1]) + Cm[i-1]
            break
    
    # return 0
    return Cm_ans

# f11
def doCx_A(w_vec, c, R, S, Omega, alpha, Lamda, theta_pitch, beta,  length, segments, k_vec ) :
    
    Cx_A = 0
    g=[]
    sum1=0
    sum2=0
    p=1
    for i in range(segments) :
        w = np.linalg.norm(w_vec[i])
        # print(w)
        k = k_vec[i]
        c1 = (w**2)*c*R / (((R*Omega)**2) * S)
        Cl, Cd = doClCd(alpha[i])
        Cl = 0.1
        Cd = 0.01
        c2 = (-Cl * sin(alpha[i]) + k*Cd*cos(alpha[i])) * sin(Lamda)  + (theta_pitch *sin(Lamda) - beta*cos(Lamda)) * (Cl*cos(alpha[i]) + k*Cd*sin(alpha[i]))
        # print(Cl,Cd)
        f1 = p*c1 * c2 / R
        # print(f1) # comment
        Cx_A += (f1*length/segments)
        p=p*-1
        #print(c1,c2) # comment
        #print(alpha[2]) # comment
    #print("Cx_A",Cx_A) # comment
    #return Cx_A
    return 0
    
    # return Cx_A

    #     g.append(f1)

    # for i in range(3,segments-3,2):
    #     sum1=sum1 + g[i]

    # for i in range(2,segments-2,2):
    #     sum2=sum2 + g[i]


    # Cx_A= length/segments/3*(g[0]+2*sum1+4*sum2+g[segments-1])

    # return Cx_A

# f12
def doCy_A(w_vec, c, R, S, Omega, alpha, Lamda, theta_pitch, beta,  length, segments, k_vec ) :
    Cy_A = 0
    g=[]
    sum1=0
    sum2=0
    for i in range(segments) :
        w = np.linalg.norm(w_vec[i])
        k = k_vec[i]
        c1 = (w**2)*c*R/ (((R*Omega)**2) * S)
        Cl, Cd = doClCd(alpha[i])
        Cl = 0.1
        Cd = 0.01
        c2 = -(-Cl * sin(alpha[i]) + k*Cd*cos(alpha[i])) * cos(Lamda)  - (theta_pitch *cos(Lamda) + beta*sin(Lamda)) * (Cl*cos(alpha[i]) + k*Cd*sin(alpha[i]))
        f1 = c1 * c2 / R
        Cy_A += f1*length/segments
        #print(c1,c2)
    # print("Cy_A",Cy_A)
    #return 0 # comment
    return Cy_A

    #     g.append(f1)

    # for i in range(3,segments-3,2):
    #     sum1=sum1 + g[i]

    # for i in range(2,segments-2,2):
    #     sum2=sum2 + g[i]


    # Cy_A= length/segments/3*(g[0]+2*sum1+4*sum2+g[segments-1])

    # return Cy_A
# f13
def doCz_A(w_vec, c, R, S, Omega, alpha, Lamda, theta_pitch, beta,  length, segments, k_vec ) :
    Cz_A = 0
    g=[]
    sum1=0
    sum2=0
    for i in range(segments) :
        w = np.linalg.norm(w_vec[i])

        k = k_vec[i]
        c1 = (w**2)*c*R/ (((R*Omega)**2) * S)
        Cl, Cd = doClCd(alpha[i])
        Cl=0.1
        Cd=0.01
        c2 = -theta_pitch*(-Cl * sin(alpha[i]) + k*Cd*cos(alpha[i])) + (Cl*cos(alpha[i]) + k*Cd*sin(alpha[i]))
        f1 = c1 * c2 / R
        Cz_A += f1*length/segments
    # print("Cz_A",Cz_A)
    return 0 # comment
    # return Cz_A

    #     g.append(f1)

    # for i in range(3,segments-3,2):
    #     sum1=sum1 + g[i]

    # for i in range(2,segments-2,2):
    #     sum2=sum2 + g[i]


    # Cz_A= length/segments/3*(g[0]+2*sum1+4*sum2+g[segments-1])

    # return Cz_A
# f8
def doCm_x_A(w_vec, c, R, S, Omega, alpha, Lamda, theta_pitch, beta, length, segments, k_vec) :
    Cm_x_A = 0
    g=[]
    sum1=0
    sum2=0
    for i in range(segments) :
        w = np.linalg.norm(w_vec[i])
        k = k_vec[i]
        c1 = (w**2)*c*R/ (((R*Omega)**2) * S)
        Cl, Cd = doClCd(alpha[i])
        Cm = doCm(alpha[i])
        eta = (i+0.5)*length / segments
        c2 = (Cl*sin(alpha[i]) + k*Cd*cos(alpha[i]))*sin(Lamda)*eta/ R - (theta_pitch*sin(Lamda) - beta*cos(Lamda)) *(-Cl*cos(alpha[i]) + k*Cd*sin(alpha[i]))*eta / R + (c/R)*Cm*cos(Lamda)
        f1 = c1  * c2 / R
        Cm_x_A += f1*length/segments
    # return 0
    # return Cm_x_A
    #     g.append(f1)
    
    # for i in range(3,segments-3,2):
    #     sum1=sum1 + g[i]

    # for i in range(2,segments-2,2):
    #     sum2=sum2 + g[i]


    # Cm_x_A= length/segments/3*(g[0]+2*sum1+4*sum2+g[segments-1])
    return 0 # comment
    # return Cm_x_A
# f9
def doCm_y_A(w_vec, c, R, S, Omega, alpha, Lamda, theta_pitch, beta,  length, segments, x_ac, k_vec) :
    integral_term = 0
    g=[]
    sum1=0
    sum2=0

    for i in range(segments) :
        w = np.linalg.norm(w_vec[i])
        k = k_vec[i]
        c1 = (w**2)*c*R/ (((R*Omega)**2) * S)
        Cl, Cd = doClCd(alpha[i])
        Cm = doCm(alpha[i])
        eta = (i+0.5)*length / segments
        c2 = -(Cl*sin(alpha[i]) + k*Cd*cos(alpha[i]))*cos(Lamda)*eta/ R + (theta_pitch*cos(Lamda) + beta*sin(Lamda)) *(-Cl*cos(alpha[i]) + k*Cd*sin(alpha[i]))*eta / R - (c/R)*Cm*cos(Lamda)
        f1 = c1 * c2 / R
      
        integral_term += f1*length/segments

    Cz = doCz_A(w_vec, c, R, S, Omega, alpha, Lamda, theta_pitch, beta,  length, segments, k_vec )
    sum_term = x_ac / R * Cz
    return 0
    # return integral_term - sum_term

    #     Cm_y_A += f1*length/segments

    # return Cm_y_A
    #     g.append(f1)

    # for i in range(3,segments-3,2):
    #     sum1=sum1 + g[i]

    # for i in range(2,segments-2,2):
    #     sum2=sum2 + g[i]


    # integral_term= length/segments/3*(g[0]+2*sum1+4*sum2+g[segments-1])

    # Cz = doCz_A(w_vec, c, R, S, Omega, alpha, Lamda, theta_pitch, beta,  length, segments, k_vec )
    # sum_term = x_ac / R * Cz
    # return integral_term - sum_term

# f10
def doCm_z_A(w_vec, c, R, S, Omega, alpha, Lamda, theta_pitch, beta,  length, segments,x_ac, k_vec) :
    integral_term = 0
    g=[]
    sum1=0
    sum2=0

    for i in range(segments) :
        w = np.linalg.norm(w_vec[i])
        k = k_vec[i]
        c1 = (w**2)*c*R/ (((R*Omega)**2) * S)
        Cl, Cd = doClCd(alpha[i])
        Cm = doCm(alpha[i])
        eta = (i+0.5)*length / segments
        c2 = -theta_pitch*(Cl*sin(alpha[i]) + k*Cd*cos(alpha[i]))*eta/ R - (-Cl*cos(alpha[i]) + k*Cd*sin(alpha[i]))*eta / R + (c/R)*Cm*cos(Lamda)
        f1 = c1 * c2 / R

        integral_term += f1*length/segments

    Cy = doCy_A(w_vec, c, R, S, Omega, alpha, Lamda, theta_pitch, beta,  length, segments, k_vec )

    sum_term = x_ac / R * Cy

    return 0
    # return integral_term + sum_term

    #     g.append(f1)

    # for i in range(3,segments-3,2):
    #     sum1=sum1 + g[i]

    # for i in range(2,segments-2,2):
    #     sum2=sum2 + g[i]

    # integral_term= length/segments/3*(g[0]+2*sum1+4*sum2+g[segments-1])

    # Cy = doCy_A(w_vec, c, R, S, Omega, alpha, Lamda, theta_pitch, beta,  length, segments, k_vec )

    # sum_term = x_ac / R * Cy

    # print(g[segments-1])
    # return g[segments-1]

    # return integral_term + sum_term

# f2
def doCx_G(phi, theta, psi, Phi0, Theta0, Psi0) :

    return (-cos(theta)*cos(psi)*sin(Theta0) + cos(theta)*sin(psi)*sin(Phi0)*cos(Theta0) - sin(theta)*cos(Phi0)*cos(Theta0))
# f3
def doCy_G(phi, theta, psi, Phi0, Theta0, Psi0) :

    return ((cos(phi)*sin(psi)-sin(phi)*sin(theta)*cos(psi)) * sin(Theta0) + (cos(phi)*cos(psi)+sin(phi)*sin(theta)*sin(psi))*sin(Phi0)*cos(Theta0) + sin(phi)*cos(theta)*cos(Phi0)*cos(Theta0))
# f4
def doCz_G(phi, theta, psi, Phi0, Theta0, Psi0) :

    return -(sin(phi)*sin(psi)+cos(phi)*sin(theta)*cos(psi)) * sin(Theta0) - (sin(phi)*cos(psi)-cos(phi)*sin(theta)*sin(psi))*sin(Phi0)*cos(Theta0) + cos(phi)*cos(theta)*cos(Phi0)*cos(Theta0)



    # Cx_A = coefficients.doCx_A(w_vec_1, c, R, S, Omega, alpha_vec_1, Lamda_1, theta_1, beta_1,  length_1, segments, k_vec_1)
    # # # print("Cx_A",Cx_A) ##
    # Cx_A += coefficients.doCx_A(w_vec_2, c, R, S, Omega, alpha_vec_2, Lamda_2, theta_2, beta_2,  length_2, segments, k_vec_2)
    # # print("Cx_A",Cx_A) ##

    # Cy_A = coefficients.doCy_A(w_vec_1, c, R, S, Omega, alpha_vec_1, Lamda_1, theta_1, beta_1,  length_1, segments, k_vec_1)
    # # # print("Cy_A",Cy_A)
    # Cy_A += coefficients.doCy_A(w_vec_2, c, R, S, Omega, alpha_vec_2, Lamda_2, theta_2, beta_2,  length_2, segments, k_vec_2)
    # #print("Cy_A",Cy_A) ##

    # Cz_A = coefficients.doCz_A(w_vec_1, c, R, S, Omega, alpha_vec_1, Lamda_1, theta_1, beta_1,  length_1, segments, k_vec_1)
    # # # #print("Cz_A",Cz_A)
    # Cz_A += coefficients.doCz_A(w_vec_2, c, R, S, Omega, alpha_vec_2, Lamda_2, theta_2, beta_2,  length_2, segments, k_vec_2)
    # # # #print(coefficients.doCz_A(w_vec_2, c, R, S, Omega, alpha_vec_2, Lamda_2, theta_2, beta_2,  length_2, segments))
    # #print("Cz_A",Cz_A) ##

    # Cx_G = coefficients.doCx_G(phi, theta, psi, Phi0, Theta0, Psi0)
    # print("Cx_g",Cx_G)
    # Cy_G = coefficients.doCy_G(phi, theta, psi, Phi0, Theta0, Psi0)
    # print("Cy_g",Cy_G)
    # Cz_G = coefficients.doCz_G(phi, theta, psi, Phi0, Theta0, Psi0)
    # print("Cz_g",Cz_G)

    # # Moment coefficients
    # Cm_x_A = coefficients.doCm_x_A(w_vec_1, c, R, S, Omega, alpha_vec_1, Lamda_1, theta_1, beta_1,  length_1, segments, k_vec_1)
    # # # #print("Cm_x_A",Cm_x_A)
    # Cm_x_A += coefficients.doCm_x_A(w_vec_2, c, R, S, Omega, alpha_vec_2, Lamda_2, theta_2, beta_2,  length_2, segments, k_vec_1)
    # #print("Cm_x_A",Cm_x_A)
    # Cm_y_A = coefficients.doCm_y_A(w_vec_1, c, R, S, Omega, alpha_vec_1, Lamda_1, theta_1, beta_1,  length_1, segments, x_ac_1, k_vec_1)
    # # # #print("Cm_y_A",Cm_y_A)
    # Cm_y_A += coefficients.doCm_y_A(w_vec_2, c, R, S, Omega, alpha_vec_2, Lamda_2, theta_2, beta_2,  length_2, segments, x_ac_2, k_vec_1)
    # #print("Cm_y_A",Cm_y_A)
    # Cm_z_A = coefficients.doCm_z_A(w_vec_1, c, R, S, Omega, alpha_vec_1, Lamda_1, theta_1, beta_1,  length_1, segments, x_ac_1, k_vec_1)
    # # # #print("Cm_z_A",Cm_z_A)
    # Cm_z_A += coefficients.doCm_z_A(w_vec_2, c, R, S, Omega, alpha_vec_2, Lamda_2, theta_2, beta_2,  length_2, segments, x_ac_2, k_vec_1)
    # #print("Cm_z_A",Cm_z_A)

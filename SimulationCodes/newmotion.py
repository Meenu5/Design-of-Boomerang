import numpy as np
from numpy.linalg import inv
from numpy.linalg import solve
import pandas as pd
import time
import copy
import coefficients
import kinematics
import transformations
from scipy.integrate import RK45

tan = np.tan
cos = np.cos
sin = np.sin
pi = np.pi
transpose = np.transpose

### PARAMETERS ###
c = 0.0488 # chord
R = 0.269 # Radius
S = 0.228 # Surface Area
m = 0.13 # mass
rho = 1.225
g = 9.8

length_1 = 0.3
length_2 = 0.3
theta_1 = 0.
theta_2 = 0.
Lamda_1 = pi*2/3
Lamda_2 = pi*4/3
x_ac_1 = 0.0723
x_ac_2 = 0.0723
beta_1 = 0.
beta_2 = 0.
segments = 30

inertia_jframe = np.array([[1.9*10**-3, 0., 0.],[0,4.88*10**-6, 0.],[0.,0.,1.9*10**-3]])
x_ac_vec = np.array([x_ac_1,x_ac_2])
length_vec = np.array([length_1, length_2]) # length of blades
theta_pitch_vec = np.array([theta_1, theta_2]) # theta_pitch of blades
Lamda_vec = np.array([Lamda_1,Lamda_2]) # Lamda of blades
beta_vec = np.array([0.,0.]) # beta of blades

delta_t = 0.0001 # time step
total_time = 0.01
total_steps = total_time/delta_t
curr_time = 0

# Initial Conditions
u_vec_d = np.array([0.,0.,0.])
u_vec = np.array([10.,0.,0.]) # initialize u_vec
omega_vec = np.array([0.,0.,60]) # initialize omega_vec
Omega = copy.copy(omega_vec[2]) # r at t = 0

path_vec = np.array([0.,0.,1.]) # initial X0,Y0,Z0
phi = 60.*pi/180.
theta = 0.
psi = 0.

Phi0 = 60.*pi/180.
Theta0 = 0.
Psi0 = 0.

Phi = 60.*pi/180.
Psi = 0.
Theta = 0.
lamda = 0.
### END PARAMETERS ####

### OVERALL TRAJECTORY ####
trajectory_u_vec = []
trajectory_omega_vec = []
trajectory_Phi_Theta_Psi = []
trajectory_path = []
trajectory_phi_theta_psi = []


# Unchanged params
Tj1 = transformations.doT0TjTransformation(Lamda_vec[0],theta_pitch_vec[0],beta_vec[0]) # blade 1
Tj2 = transformations.doT0TjTransformation(Lamda_vec[1],theta_pitch_vec[1],beta_vec[1]) # blade 2

r_j_vec_1 = kinematics.doPositionVector(Tj1,length_vec[0],x_ac_vec[0],segments)
r_j_vec_2 = kinematics.doPositionVector(Tj2,length_vec[1],x_ac_vec[1],segments)

J1 = np.matmul(transpose(Tj1),inertia_jframe)
J2 = np.matmul(transpose(Tj2),inertia_jframe)
J = J1+J2


# for i in range(3) :
#     for j in range(3) :
#         if J[i][j] < 10**-5 :
#             J[i][j] = 0


# initial velocity
trajectory_u_vec.append(copy.copy(u_vec))
trajectory_omega_vec.append(copy.copy(omega_vec))
trajectory_Phi_Theta_Psi.append(copy.copy(np.array([Phi, Theta, Psi])))
trajectory_phi_theta_psi.append(copy.copy(np.array([Phi,Theta,Psi])))
trajectory_path.append(copy.copy(path_vec))


a = time.time()
steps = 0
while steps <= total_steps :
    ## Changes over time ###

    # kinematics
    v_vec_1 = kinematics.doVelBlade(u_vec, omega_vec, r_j_vec_1)
    v_vec_2 = kinematics.doVelBlade(u_vec, omega_vec, r_j_vec_2)

    w_vec_1 = kinematics.doRelativeAirVelBlade(v_vec_1,Tj1)
    w_vec_2 = kinematics.doRelativeAirVelBlade(v_vec_2,Tj2)

   
    alpha_vec_1, k_vec_1 = kinematics.doAlpha(w_vec_1)
    alpha_vec_2, k_vec_2 = kinematics.doAlpha(w_vec_2)





    # Coefficients - Aerodynamic and Body
    Cx_A = coefficients.doCx_A(w_vec_1, c, R, S, Omega, alpha_vec_1, Lamda_1, theta_1, beta_1,  length_1, segments, k_vec_1)
    # print("Cx_A",Cx_A) ##
    Cx_A += coefficients.doCx_A(w_vec_2, c, R, S, Omega, alpha_vec_2, Lamda_2, theta_2, beta_2,  length_2, segments, k_vec_2)
    # print("Cx_A",Cx_A) ##

    Cy_A = coefficients.doCy_A(w_vec_1, c, R, S, Omega, alpha_vec_1, Lamda_1, theta_1, beta_1,  length_1, segments, k_vec_1)
    # print("Cy_A",Cy_A)
    Cy_A += coefficients.doCy_A(w_vec_2, c, R, S, Omega, alpha_vec_2, Lamda_2, theta_2, beta_2,  length_2, segments, k_vec_2)
    # print("Cy_A",Cy_A) ##

    Cz_A = coefficients.doCz_A(w_vec_1, c, R, S, Omega, alpha_vec_1, Lamda_1, theta_1, beta_1,  length_1, segments, k_vec_1)
    # # #print("Cz_A",Cz_A)
    Cz_A += coefficients.doCz_A(w_vec_2, c, R, S, Omega, alpha_vec_2, Lamda_2, theta_2, beta_2,  length_2, segments, k_vec_2)
    # # #print(coefficients.doCz_A(w_vec_2, c, R, S, Omega, alpha_vec_2, Lamda_2, theta_2, beta_2,  length_2, segments))
    #print("Cz_A",Cz_A) ##

    Cx_G = coefficients.doCx_G(phi, theta, psi, Phi0, Theta0, Psi0)
    # print("Cx_g",Cx_G) # comment  
    Cy_G = coefficients.doCy_G(phi, theta, psi, Phi0, Theta0, Psi0)
    # print("Cy_g",Cy_G) # comment
    Cz_G = coefficients.doCz_G(phi, theta, psi, Phi0, Theta0, Psi0)
    # print("Cz_g",Cz_G) # comment

    # Moment coefficients
    Cm_x_A = coefficients.doCm_x_A(w_vec_1, c, R, S, Omega, alpha_vec_1, Lamda_1, theta_1, beta_1,  length_1, segments, k_vec_1)
    # # #print("Cm_x_A",Cm_x_A)
    Cm_x_A += coefficients.doCm_x_A(w_vec_2, c, R, S, Omega, alpha_vec_2, Lamda_2, theta_2, beta_2,  length_2, segments, k_vec_1)
    #print("Cm_x_A",Cm_x_A)
    Cm_y_A = coefficients.doCm_y_A(w_vec_1, c, R, S, Omega, alpha_vec_1, Lamda_1, theta_1, beta_1,  length_1, segments, x_ac_1, k_vec_1)
    # # #print("Cm_y_A",Cm_y_A)
    Cm_y_A += coefficients.doCm_y_A(w_vec_2, c, R, S, Omega, alpha_vec_2, Lamda_2, theta_2, beta_2,  length_2, segments, x_ac_2, k_vec_1)
    #print("Cm_y_A",Cm_y_A)
    Cm_z_A = coefficients.doCm_z_A(w_vec_1, c, R, S, Omega, alpha_vec_1, Lamda_1, theta_1, beta_1,  length_1, segments, x_ac_1, k_vec_1)
    # # #print("Cm_z_A",Cm_z_A)
    Cm_z_A += coefficients.doCm_z_A(w_vec_2, c, R, S, Omega, alpha_vec_2, Lamda_2, theta_2, beta_2,  length_2, segments, x_ac_2, k_vec_1)
    #print("Cm_z_A",Cm_z_A)

    # # Linear Acceleration Calculation
    u_vec_d[0]  = 0.5*rho*((R*Omega)**2)*S*Cx_A / m- g*Cx_G - omega_vec[1]*u_vec[2] + omega_vec[2]*u_vec[1]
    u_vec_d[1]  = 0.5*rho*((R*Omega)**2)*S*Cy_A / m- g*Cy_G - omega_vec[2]*u_vec[0] + omega_vec[0]*u_vec[2]
    u_vec_d[2]  = 0.5*rho*((R*Omega)**2)*S*Cz_A / m- g*Cz_G - omega_vec[0]*u_vec[1] + omega_vec[1]*u_vec[0]


    # # Rate of angular velocity calculation

    # t1 = (J[2][2]-J[1][1])*omega_vec[1]*omega_vec[2]
    # t2 = (J[0][1]*omega_vec[0] + J[1][2]*omega_vec[2])*omega_vec[2]
    # t3 = -(J[0][2]*omega_vec[0] + J[1][2]*omega_vec[1])*omega_vec[1]

    # rhs1 = 0.5*rho*((R*Omega)**2)*R*S*Cm_x_A - t1 - t2 - t3
    # lhs1 = np.array([J[0][0],-J[0][1],-J[0][2]])


    # t1 = (J[0][0]-J[2][2])*omega_vec[0]*omega_vec[2]
    # t2 = (J[1][2]*omega_vec[1] + J[0][2]*omega_vec[0])*omega_vec[0]
    # t3 = -(J[0][1]*omega_vec[0] + J[0][2]*omega_vec[2])*omega_vec[2]

    # rhs2 = 0.5*rho*((R*Omega)**2)*R*S*Cm_y_A - t1 - t2 - t3
    # lhs2 = np.array([J[1][1],-J[1][2],-J[0][1]])


    # t1 = (J[1][1]-J[0][0])*omega_vec[0]*omega_vec[1]
    # t2 = (J[0][2]*omega_vec[2] + J[0][1]*omega_vec[1])*omega_vec[1]
    # t3 = -(J[1][2]*omega_vec[2] + J[0][1]*omega_vec[0])*omega_vec[0]

    # rhs3 = 0.5*rho*((R*Omega)**2)*R*S*Cm_z_A - t1 - t2 - t3
    # lhs3 = np.array([J[2][2],-J[0][2],-J[1][2]])

    # lhs = np.array([lhs1, lhs2, lhs3])
    # rhs = np.array([rhs1, rhs2, rhs3])
    # omega_vec_d = np.linalg.solve(lhs,rhs)

    # Linear Acceleration Calculation
    c1 = 0.5*rho*((R*Omega)**2)*S*Cx_A / m- g*Cx_G
    c2 = 0.5*rho*((R*Omega)**2)*S*Cy_A / m- g*Cy_G
    c3 = 0.5*rho*((R*Omega)**2)*S*Cz_A / m- g*Cz_G
    c4 = 0.5*rho*((R*Omega)**2)*R*S*Cm_x_A
    c5 = 0.5*rho*((R*Omega)**2)*R*S*Cm_y_A
    c6 = 0.5*rho*((R*Omega)**2)*R*S*Cm_z_A

    # rk45
    def funvel(t,x) :
        global J, c1,c2,c3,c4,c5,c6 
        
        x0d  = c1 - x[4]*x[2] + x[5]*x[1]
        x1d  = c2 - x[5]*x[0] + x[3]*x[2]
        x2d  = c3 - x[3]*x[1] + x[4]*x[0]

        x3d  = (-J[2][2]*x[3]*x[5]-J[0][1]*x[3]*x[5]-c5)/J[0][1]
        x4d  = (J[2][2]*x[4]*x[5]+J[0][1]*x[3]*x[5]-c4)/J[0][1]
        x5d  = (c6-J[0][1]*x[4]*x[4]-J[0][1]*x[3]*x[3])/J[2][2]

        return np.array([x0d,x1d,x2d,x3d,x4d,x5d])

    obj1 = RK45(funvel,curr_time,np.array([u_vec[0],u_vec[1],u_vec[2],omega_vec[0],omega_vec[1],omega_vec[2]]),curr_time+delta_t)
    obj1.step()
    obj2 = obj1.dense_output()
    ans_vel = obj2.__call__(curr_time+delta_t)

    # computing lamda and r_n angular rate of Z axis of nonspinning frame
    ### Uncomment
    r_n = omega_vec[2] + (u_vec_d[1]*u_vec[0] - u_vec_d[0]*u_vec[1]) / (u_vec[0]**2 + u_vec[1]**2)
    lamda += (omega_vec[2] - r_n)*delta_t

    def funcapangle(t,x) :
        global omega_vec,lamda,r_n
        Theta = x[1]
        Phi = x[2]
        Psi_d = (omega_vec[0]*sin(lamda) + omega_vec[1]*cos(lamda)) * sin(Phi)/cos(Theta) + r_n*cos(Phi)/cos(Theta)
        Theta_d = (omega_vec[0]*sin(lamda) + omega_vec[1]*cos(lamda))*cos(Phi) - r_n*sin(Phi)
        Phi_d = (omega_vec[0]*cos(lamda)-omega_vec[1]*sin(lamda)) + (omega_vec[0]*sin(lamda)+omega_vec[1]*cos(lamda))*sin(Phi)*tan(Theta) + r_n*cos(Phi)*tan(Theta)

        return np.array([Psi_d,Theta_d,Phi_d])

    obj1 = RK45(funcapangle,curr_time,np.array([Psi,Theta,Phi]),curr_time+delta_t)
    obj1.step()
    obj2 = obj1.dense_output()
    ans_capangles = obj2.__call__(curr_time+delta_t)

    def funangle(t,x) :
        global omega_vec,lamda,r_n
        theta = x[1]
        phi = x[2]
        psi_d = omega_vec[1]*sin(phi) / cos(theta) + omega_vec[2] * cos(phi) / cos(theta)
        theta_d = omega_vec[1] * cos(phi) - omega_vec[2]*sin(phi)
        phi_d = omega_vec[0] + omega_vec[1]*sin(phi)* tan(theta) + omega_vec[2]*cos(phi)*tan(theta)

        return np.array([psi_d,theta_d,phi_d])

    obj1 = RK45(funangle,curr_time,np.array([psi,theta,phi]),curr_time+delta_t)
    obj1.step()
    obj2 = obj1.dense_output()
    ans_angles = obj2.__call__(curr_time+delta_t)

    
    # next time step

    u_vec = ans_vel[:3]
    omega_vec = ans_vel[3:]

    Psi = ans_capangles[0]
    Theta = ans_capangles[1]
    Phi = ans_capangles[2]

    psi = ans_angles[0]
    theta = ans_angles[1]
    phi = ans_angles[2]

    T0 = transformations.doT0Transformation(phi, theta, psi)
    Ti = transformations.doTiTransformation(Phi0, Theta0, Psi0)
    U = np.matmul(np.matmul(inv(Ti),inv(T0)),u_vec)

    path_vec += U*delta_t

    Angle_vec = np.array([Phi, Theta, Psi])
    angle_vec = np.array([phi, theta, psi])

    trajectory_u_vec.append(copy.copy(u_vec))
    trajectory_omega_vec.append(copy.copy(omega_vec))
    trajectory_Phi_Theta_Psi.append(copy.copy(Angle_vec))
    trajectory_phi_theta_psi.append(copy.copy(angle_vec))
    trajectory_path.append(copy.copy(path_vec))

    # #print(alpha_vec_1)
    # #print(alpha_vec_2)
    print("Step ",steps)
    steps += 1
    curr_time += delta_t

 

    if steps%1 == 0 :
        df_u_vec = pd.DataFrame(trajectory_u_vec, dtype=None, copy=False)
        df_omega_vec = pd.DataFrame(trajectory_omega_vec, dtype=None, copy=False)
        df_Phi_Theta_Psi = pd.DataFrame(trajectory_Phi_Theta_Psi, dtype=None, copy=False)
        df_phi_theta_psi = pd.DataFrame(trajectory_phi_theta_psi, dtype=None, copy=False)
        df_path = pd.DataFrame(trajectory_path, dtype=None, copy=False)

        df_u_vec.to_csv("linearVel.csv")
        df_omega_vec.to_csv("angularVel.csv")
        df_Phi_Theta_Psi.to_csv("eulerAnglesCap.csv")
        df_phi_theta_psi.to_csv("eulerAngles.csv")
        df_path.to_csv("coordinates.csv")

b = time.time()

print("Time taken in seconds",b-a)
# print(alpha_vec_1)
# print(alpha_vec_2)
# # print(trajectory_u_vec)
# # print(trajectory_omega_vec)
# # print(trajectory_Phi_Theta_Psi)
# # print(trajectory_phi_theta_psi)
# # print(trajectory_path)

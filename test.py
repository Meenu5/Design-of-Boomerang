import coefficients
import numpy as np 

Cl,Cd = coefficients.doClCd(0)




def doCx_A(w_vec, c, R, S, Omega, alpha, Lamda, theta_pitch, beta,  length, segments ) :
    Cx_A = 0
    for i in range(segments) :
        w = np.linalg.norm(w_vec[i])
        c1 = (w**2)*c*R / (((R*Omega)**2) * S)
        Cl, Cd = doClCd(alpha[i])
        print("Cl",Cl)
        c2 = -Cl * sin(alpha[i]) * sin(Lamda)  + (theta_pitch *sin(Lamda) - beta*cos(Lamda)) * (Cl*cos(alpha[i]))
        f1 = c1 * c2 / R
        Cx_A += f1*length/segments

    return Cx_A

print(Cx_A)



print("first",-theta_pitch*Cl*sin(alpha[i])*eta/ R)
print("last but one",alpha[i],Cl,eta,R,- (-Cl*cos(alpha[i])*eta / R))
print("last",(c/R)*Cm*cos(Lamda))
print("c1,c2",c1,c2 )

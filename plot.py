import coefficients
import numpy as np 
import matplotlib.pyplot as plt 
pi = np.pi

alpha = -pi
alpha_array = []
Cl_array = []
Cd_array = []
Cm_array = []

# while(alpha<=pi) :
#     Cl,Cd = coefficients.doClCd(alpha)
#     Cl_array += [Cl]
#     Cd_array += [Cd]
#     Cm = coefficients.doCm(alpha)
#     Cm_array += [Cm]
#     alpha_array += [alpha]
#     alpha += 0.05
    
Cl,Cd = coefficients.doClCd(0.20944)

print(Cl,'reoblern')
plt.plot(alpha_array,Cm_array)
plt.xlabel("Angle of attack in radians")
plt.ylabel("Moment Coefficient Cm")
plt.title("Cm vs alpha")
plt.show()
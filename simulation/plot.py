import coefficients
import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd 
pi = np.pi

# alpha = -pi
# alpha_array = []
# Cl_array = []
# Cd_array = []
# Cm_array = []

# # while(alpha<=pi) :
# #     Cl,Cd = coefficients.doClCd(alpha)
# #     Cl_array += [Cl]
# #     Cd_array += [Cd]
# #     Cm = coefficients.doCm(alpha)
# #     Cm_array += [Cm]
# #     alpha_array += [alpha]
# #     alpha += 0.05
    
# Cl,Cd = coefficients.doClCd(0.20944)

# print(Cl)
# plt.plot(alpha_array,Cm_array)
# plt.xlabel("Angle of attack in radians")
# plt.ylabel("Moment Coefficient Cm")
# plt.title("Cm vs alpha")
# plt.show()


df = pd.read_csv("coordinates.csv")
x = list(df["0"])
y = list(df["2"])
y_the = [-i**2/200*9.81 for i in x]
# plt.plot(x,y)
x = x[1:]
y = y[1:]
ratio = []
time = []
for i in range(len(x)) :
    ratio += [(x[i]**2)/y[i]]
    time += [(i+1)*0.01]
    gt += 

plt.plot(time,y)
plt.plot(time,y_the[1:])

plt.xlabel("X coordinate")
plt.ylabel("Y coordinate")
plt.title("X vs Y")
plt.show()
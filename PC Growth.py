import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

t0 = float(0)
tf = float(1000)
dt = 1
n = int(tf / dt)

initialC = float(input("initial C"))
initialP = float(input("initial P"))
initialR = float(input("initial R"))
initialT = float(input("initial T"))

C = float(initialC)
P = float(initialP)
R = float(initialR)
T = float(initialT)

Co = float(10**6)
Po = float(10**5)
kc = float(7 * 10**2)
muc = float((20 * kc) / Po)
Kc = float(.1)
lc = float(10**-7)
kp = float(.2)
mup = float(20 * kp)
Kp = float(Co / 100)
lp = float(.15)
kr = float(.2)
lr = float(.22)
nup = float((.02 * lr) / (Po * (1 - lp / kp)))
nuc = nup
kt = float(3300)
Kt = Kc
lt = float(.3)

x = np.linspace (t0 , tf , tf / dt)
C = np.zeros([n])
C[0] = initialC
for i in range(1, n):
    C[i] = dt * [(kc + muc * P) * C * (1 - (C / Co)) - (lc * C * T) / (Kc + (1 - R))] + C[float(i) - 1]
plt.plot(x, C, 'o-')
    
P = np.zeros([n])
P[0] = initialP
for i in range(1, n):
    P[i] = dt * [(kp + (mup * C) / (Kp + C)) * P * (1 - (P / Po)) - (lp * P)] + P[float(i) - 1]
plt.plot(x, P, 'o-')
    
R = np.zeros([n])
R[0] = initialR
for i in range(1, n):
    R[i] = dt * [kr - (lr + (nup * P) + (nuc * C)) * R] + R[float(i) - 1]
plt.plot(x, R, 'o-')
      
T = np.zeros([n])
T[0] = initialT
for i in range(1, n):
    T[i] = dt * [((kt * R) / (Kt + (1 - R))) - (lt * T)] + T[float(i) - 1]
plt.plot(x, T, 'o-')
    
plt.show()
    

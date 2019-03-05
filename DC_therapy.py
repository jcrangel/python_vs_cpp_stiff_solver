from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp 
import time

# PARAM
aH = 1e-4
MuH = 0.005
rH = 10e-2
KH = 1
aC = 1e-4
MuC = 0.01925
rC = 0.00004e-2
KC = 1
KT = 1e12
MuD = 0.009625
rI = 1e-2
MuIC = 1e-7
MuI = 1e-2
rT = 0.002
aT = 0.1136
eT = 50
hT = 5.2e5
aTBeta = 0.69
eTBeta = 1e4
rTBeta = 5.57e-6
MuBeta = 6.93
aGammaC = 1.02e-4
MuGamma = 0.102
gMl = 1.44
aMlGamma = 2.89
eMlGamma = 3.38e5
MuMl = 0.0144


def DC_model(t,y):
    T,H,CTL,Den,IL2,FBeta,FGamma,Ml = y

    dT = rT * T * np.log(KT / T) - (aT * T * CTL * Ml / (eT + Ml)) * ((aTBeta * FBeta + eTBeta) / (FBeta + eTBeta))
    dH = aH - MuH * H + rH * Den*(H*(1 - H / KH))
    dC = aC - MuC * CTL + rC * IL2*(CTL*(1 - CTL / KC))
    dDen = -MuD * Den * CTL
    dI = rI * H * Den - MuIC * CTL * IL2 - MuI * IL2
    dFbeta = rTBeta * T - MuBeta * FBeta
    dFgamma = aGammaC * CTL - MuGamma * FGamma
    dMl = gMl + (aMlGamma * FGamma) / (FGamma + eMlGamma) - MuMl * Ml

    return np.array([dT,dH,dC,dDen,dI,dFbeta,dFgamma,dMl])

def DC_model_jacobian(t,y):
    T,H,CTL,Den,IL2,FBeta,FGamma,Ml = y
    M = np.empty((8,8))
    M[0, 0] = -rT + rT * np.log(KT / T) - (aT*CTL*(eTBeta + aTBeta * FBeta) * Ml) / ((eTBeta + FBeta)* (eT + Ml))
    M[0, 1] = 0
    M[0, 2] = -(aT * (eTBeta + aTBeta * FBeta)* Ml * T) / ((eTBeta + FBeta) * (eT + Ml))
    M[0, 3] = 0
    M[0, 4] = 0
    M[0, 5] = -((aT*aTBeta*CTL*Ml*T) / ((eTBeta + FBeta) *(eT + Ml))) + (aT*CTL *(eTBeta + aTBeta * FBeta)* Ml*T) / (np.power(eTBeta + FBeta, 2) * (eT + Ml))
    M[0, 6] = 0
    M[0, 7] = (aT * CTL*(eTBeta + aTBeta * FBeta)* Ml * T) / ((eTBeta + FBeta)* np.power(eT + Ml, 2)) - (aT * CTL*(eTBeta + aTBeta * FBeta)* T) / ((eTBeta + FBeta)*(eT + Ml))
    #row 2
    M[1, 0] = 0
    M[1, 1] = -MuH - (rH * Den * H) / KH + rH * Den* (1 - H / KH)
    M[1, 2] = 0
    M[1, 3] = rH * H*(1 - H / KH)
    M[1, 4] = 0
    M[1, 5] = 0
    M[1, 6] = 0
    M[1, 7] = 0
    #row 3
    M[2, 0] = 0
    M[2, 1] = 0
    M[2, 2] = -MuC - (rC * CTL * IL2) / KC + rC * (1 - CTL / KC)*IL2
    M[2, 3] = 0
    M[2, 4] = rC * CTL*(1 - CTL / KC)
    M[2, 5] = 0
    M[2, 6] = 0
    M[2, 7] = 0
    #row 4
    M[3, 0] = 0
    M[3, 1] = 0
    M[3, 2] = -MuD * Den
    M[3, 3] = -MuD * CTL
    M[3, 4] = 0
    M[3, 5] = 0
    M[3, 6] = 0
    M[3, 7] = 0
    #row 5
    M[4, 0] = 0
    M[4, 1] = rI * Den
    M[4, 2] = -MuIC * IL2
    M[4, 3] = rI * H
    M[4, 4] = -MuI - MuIC * CTL
    M[4, 5] = 0
    M[4, 6] = 0
    M[4, 7] = 0
    #row 6
    M[5, 0] = rTBeta
    M[5, 1] = 0
    M[5, 2] = 0
    M[5, 3] = 0
    M[5, 4] = 0
    M[5, 5] = -MuBeta
    M[5, 6] = 0
    M[5, 7] = 0

    #row 7
    M[6, 0] = 0
    M[6, 1] = 0
    M[6, 2] = aGammaC
    M[6, 3] = 0
    M[6, 4] = 0
    M[6, 5] = 0
    M[6, 6] = -MuGamma
    M[6, 7] = 0

    #row 8
    M[7, 0] = 0
    M[7, 1] = 0
    M[7, 2] = 0
    M[7, 3] = 0
    M[7, 4] = 0
    M[7, 5] = 0
    M[7, 6] = -((aMlGamma * FGamma) / np.power(eMlGamma + FGamma, 2)) + aMlGamma / (eMlGamma + FGamma)
    M[7, 7] = -MuMl



y0 = [6e4, 0, 0, 0, 0, 0, 0, 0]
t0 = 0.0 
delay = 232
tf = 168 + delay
ef=0.05

start = time.time()
sol = solve_ivp(DC_model,[t0,tf],y0,method='LSODA',atol=1e-4,jac=DC_model_jacobian)  #,method='LSODA'
plt.plot(sol.t, sol.y[0,:], 'b')
# plt.plot(sol.t, sol.y[:, 1], 'g', label='omega(t)')

#First injection
y0 = sol.y[:,-1] + [0, 0, 0, 1e6 * ef, 0, 0, 0, 0]
t0 = tf 
tf = tf + 168   
sol = solve_ivp(DC_model,[t0,tf],y0,method='LSODA',atol=1e-4,jac=DC_model_jacobian)  #,method='LSODA'
plt.plot(sol.t, sol.y[0,:], 'b')

#Second injection
y0 = sol.y[:,-1] + [0, 0, 0, 1e6 * ef, 0, 0, 0, 0]
t0 = tf 
tf = tf + 168   
sol = solve_ivp(DC_model,[t0,tf],y0,method='LSODA',atol=1e-4,jac=DC_model_jacobian)  #,method='LSODA'
plt.plot(sol.t, sol.y[0,:], 'b')

#Thrid injection
y0 = sol.y[:,-1] + [0, 0, 0, 1e6 * ef, 0, 0, 0, 0]
t0 = tf 
tf = 1400   
sol = solve_ivp(DC_model,[t0,tf],y0,method='LSODA',atol=1e-4,jac=DC_model_jacobian)  #,method='LSODA'
end = time.time()
print('Solving took: ' + str((end-start)/60) + ' min')
plt.plot(sol.t, sol.y[0,:], 'b', label='T(t)')



plt.legend(loc='best')
plt.xlabel('t')
plt.grid()
plt.show()

# On laptop
# LSODA took  5.526416440804799 min ,atol=1e-5, no jac
# LSODA took  5.478658803304037 min ,atol=1e-4, no jac
# LSODA took  8.489046243826548 min min ,atol=1e-4, with explicit jac!
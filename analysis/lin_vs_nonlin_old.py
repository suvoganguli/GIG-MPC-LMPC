import numpy as np
import matplotlib.pyplot as plt

def nonlin(uk, xk, T):

    xkp1 = [0.0, 0.0, 0.0, 0.0]
    xkp1[0] = xk[0] + T * xk[2] * np.sin(xk[3]) # Edot
    xkp1[1] = xk[1] + T * xk[2] * np.cos(xk[3]) # Ndot
    xkp1[2] = xk[2] + T * uk[0]                 # Vdot
    xkp1[3] = xk[3] + T * uk[1]                 # Chidot

    return xkp1


def lin(x0_k, del_u_k, del_x_k, T):

    E0_k   = x0_k[0]
    N0_k   = x0_k[1]
    V0_k   = x0_k[2]
    Chi0_k = x0_k[3]

    del_E_k   = del_x_k[0]
    del_N_k   = del_x_k[1]
    del_V_k   = del_x_k[2]
    del_Chi_k = del_x_k[3]

    del_u1_k = del_u_k[0]
    del_u2_k = del_u_k[1]

    del_E_kp1   = del_E_k + T * (np.sin(Chi0_k) * del_V_k + V0_k * np.cos(Chi0_k) * del_Chi_k)  # Edot
    del_N_kp1   = del_N_k + T * (np.cos(Chi0_k) * del_V_k - V0_k * np.sin(Chi0_k) * del_Chi_k)  # Edot
    del_V_kp1   = del_V_k + T * del_u1_k
    del_Chi_kp1 = del_Chi_k + T * del_u2_k

    x_del_kp1 = np.array([del_E_kp1, del_N_kp1, del_V_kp1, del_Chi_kp1 ])

    return x_del_kp1

# ---------------------------------------------------------
# Run simulation

T = 0.1
t = np.arange(0,10,T)

V0 = 10
chi0 = 30 * np.pi / 180

x0 = np.array([10, 20, V0, chi0])
xk_nonlin = np.zeros([4,t.size])
xk_lin = np.zeros([4,t.size])

xk_nonlin[:,0] = x0
xk_lin[:,0] = x0

u1_k = np.zeros(t.size) + 1
u2_k = np.zeros(t.size)
uk = np.array([u1_k, u2_k])

for k in range(len(t)-1):
    xkp1_nonlin     = nonlin(uk[:,k], xk_nonlin[:,k], T)
    xk_nonlin[:,k+1]  = xkp1_nonlin

    xkp1_lin     = lin(x0, uk[:,k], xk_lin[:,k], T)
    xk_lin[:,k+1]  = xkp1_lin



plt.figure(0,figsize=(6,7))

plt.subplot(411)
plt.plot(t, xk_nonlin[0,:], t, xk_lin[0,:])
plt.ylabel('E [ft]')
plt.grid('on')

plt.subplot(412)
plt.plot(t, xk_nonlin[1,:], t, xk_lin[1,:])
plt.ylabel('N [ft]')
plt.grid('on')

plt.subplot(413)
plt.plot(t, xk_nonlin[2,:], t, xk_lin[2,:])
plt.ylabel('V [fps]')
plt.grid('on')

plt.subplot(414)
plt.plot(t, xk_nonlin[3,:]*180/np.pi, t, xk_lin[3,:]*180/np.pi)
plt.ylabel('Chi [deg]')
plt.grid('on')


plt.show()








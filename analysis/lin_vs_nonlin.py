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
    del_N_kp1   = del_N_k + T * (np.cos(Chi0_k) * del_V_k - V0_k * np.sin(Chi0_k) * del_Chi_k)  # Ndot
    del_V_kp1   = del_V_k + T * del_u1_k
    del_Chi_kp1 = del_Chi_k + T * del_u2_k

    x_del_kp1 = np.array([del_E_kp1, del_N_kp1, del_V_kp1, del_Chi_kp1 ])

    return x_del_kp1


def doublet(dt,tf,t1,t2,t3,val):
    tvec = np.arange(0,tf,dt)
    yvec = 0*tvec

    idx1 = int(np.round(t1 / dt))
    idx2 = int(np.round(t2 / dt))
    idx3 = int(np.round(t3 / dt))

    yvec[idx1:idx2] = val
    yvec[idx2:idx3] = -val*0

    return tvec, yvec


def step(dt,tf,t1,val):
    tvec = np.arange(0,tf,dt)
    yvec = 0*tvec

    idx1 = int(np.round(t1 / dt))

    yvec[idx1:n] = val

    return tvec, yvec



# ---------------------------------------------------------
# Run simulation

T = 0.1
tf = 10
t = np.arange(0,tf,T)
n = len(t)

V0 = 10
chi0 = 30 * np.pi / 180
x0 = np.array([10, 20, V0, chi0])

x_nonlin = np.zeros([4,t.size])
x0_nonlin = np.zeros([4,t.size])
del_x_lin = np.zeros([4,t.size])

x_nonlin[:,0] = x0
x0_nonlin[:,0] = x0

u0 = np.zeros(2)

del_u1 = np.zeros((1,n))
#_, del_u2 = doublet(T,tf, 2.0, 3.0, 4.0, 10*np.pi/180)
_, del_u2 = step(T, tf, 2.0, 10*np.pi/180)
del_u = np.vstack([del_u1, del_u2])

u = u0[:,None]*np.zeros((2,n)) + del_u

for k in range(len(t)-1):
    xkp1_nonlin      = nonlin(u[:,k], x_nonlin[:,k], T)
    x_nonlin[:,k+1]  = xkp1_nonlin

    x0kp1            = nonlin(u0, x0_nonlin[:,k], T)
    x0_nonlin[:,k+1] = x0kp1

    del_xkp1_lin     = lin(x0, del_u[:,k], del_x_lin[:,k], T)
    #del_xkp1_lin = lin(x0_nonlin[:,k], del_u[:, k], del_x_lin[:, k], T)
    del_x_lin[:,k+1] = del_xkp1_lin


del_x_nonlin = x_nonlin - x0_nonlin


plt.figure(0,figsize=(6,7))

plt.subplot(411)
plt.plot(t, del_x_nonlin[0,:], t, del_x_lin[0,:],'--')
plt.ylabel('E [ft]')
plt.grid('on')

plt.subplot(412)
plt.plot(t, del_x_nonlin[1,:], t, del_x_lin[1,:],'--')
plt.ylabel('N [ft]')
plt.grid('on')

plt.subplot(413)
plt.plot(t, del_x_nonlin[2,:], t, del_x_lin[2,:],'--')
plt.ylabel('V [fps]')
plt.grid('on')

plt.subplot(414)
plt.plot(t, del_x_nonlin[3,:]*180/np.pi, t, del_x_lin[3,:]*180/np.pi,'--')
plt.ylabel('Chi [deg]')
plt.grid('on')


plt.figure(1)
plt.subplot(211)
plt.plot(t,u[0,:])
plt.ylabel('Vdot [fps2]')
plt.grid('on')

plt.subplot(212)
plt.plot(t,u[1,:]*180/np.pi)
plt.ylabel('Chidot [deg]')
plt.xlabel(['t [sec]'])
plt.grid('on')


plt.figure(2)
plt.plot(del_x_lin[0,:],del_x_lin[1,:])
plt.plot(del_x_nonlin[0,:],del_x_nonlin[1,:])
plt.ylabel('North [ft]')
plt.xlabel('East [ft]')
plt.grid('on')
plt.axis('equal')

plt.show()








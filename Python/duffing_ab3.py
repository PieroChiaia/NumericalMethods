from cmath import pi
from math import ceil
import numpy as np
import matplotlib.pyplot as plt

# ==============================================================================
#  DEFINE DYNAMICAL SYSTEM
def Duffing(t,x, params):
    alpha = params["alpha"]
    beta = params["beta"]
    gamma = params["gamma"]
    delta = params["delta"]
    omega = params["omega"]

    xdot = np.array([x[1], -delta*x[1]-alpha*x[0]-beta*x[0]**3+gamma*np.cos(omega*t)])

    return xdot



# ==============================================================================
#  RK4 EXPLICIT SCHEME
def RungeKutta4(f, x0, t0, tf, dt):

    t = np.arange(t0, tf, dt)
    nt = t.size

    nx = x0.size

    x=np.zeros((nx,nt))

    x[:,0] = x0

    for k in range(nt -1):
        k1 = dt*f(t[k], x[:,k])
        k2 = dt*f(t[k] + dt/2, x[:,k]+ k1/2)
        k3 = dt*f(t[k] + dt/2, x[:,k]+ k2/2)
        k4 = dt*f(t[k] + dt, x[:,k]+ k3)

        x[:,k+1] = x[:,k] + 1/6*(k1+2*k2+2*k3+k4)

    x=np.transpose(x)
    return x,t



# ==============================================================================
#  AB3 EXPLICIT SCHEME
def AdamsBashforth3(f, x0, t0, tf, dt):

    # Initializing variables
    t = np.arange(t0, tf, dt)
    nt = t.size
    nx = x0.size
    x=np.zeros((nt,nx))

    # Initial conditions and initizialization of AB3 method with Explicit Euler
    x[0,:] = x0
    x[1,:] = x[0,:] + dt * f(dt,x[0,:])
    x[2,:] = x[0,:] + dt * f(dt,x[1,:])

    for k in range(2,nt -1):
        x[k+1,:] = x[k,:] + dt/12*( 23* f(k*dt, x[k,:]) -16*f((k-1)*dt, x[k-1,:]) + 5*f((k-2)*dt, x[k-2,:]))
    
    return x,t



# ==============================================================================
#  POINCARE MAP PLOT
def PoincareMap(sol,t,omega,t_start):
    dt = t[1]-t[0]
    p_map=np.zeros((t.size,np.size(sol[0,:])))

    for k in range(len(t)):
        ind=ceil(t_start+1/dt*2*pi*k/omega)
        if ind <= len(t):
            p_map[k,:] = sol[ind,:]
        else:
            p_map[k,:] = []
            # continue
    
    # plt.plot(p_map[:,0],p_map[:,1],'bs',label="RK4")
    plt.figure(figsize=(10,6))
    title = r"$\ddot{x} +\delta \dot{x}+\alpha x +\beta x^3 = \gamma \cos(\omega t)$"
    plt.title(title)
    plt.xlim([-1.6, 1.6])
    plt.ylim([-1.25,1.25])
    plt.scatter(p_map[:,0],p_map[:,1],s=2, c="coral")
    plt.savefig('plot'+str(t_start)+'.png',dpi=300)
    plt.close()
    return



# ==============================================================================
#  NUMERICAL SOLUTION

#  Definition of Duffing Oscillator parameters and equivalent I order system
params = {"alpha": 1, "beta": 5, "gamma": 8, "delta":0.02, "omega":0.5}
f = lambda t,x : Duffing(t,x, params)

# Initial conditions and disretization
x0=np.array([0,0])
t0=0
tf=50
N=4000
dt=tf/N


# Solution by Runge-Kutta 4 - Explicit scheme
sol,t = RungeKutta4(f,x0,t0,tf,dt)
x_RK4 = sol[:,0]
xdot_RK4 = sol[:,1]


# Solution by Adams-Bashforth 3 - Explicit scheme
sol1,t = AdamsBashforth3(f,x0,t0,tf,dt)
x_AB3 = sol1[:,0]
xdot_AB3 = sol1[:,1]



# ==============================================================================
#  SEQUENCE OF POINCARE MAPS: FOR A COMPLETE ANIMATION SET
#   tf=10000   N=6000000, change constants to: 
#   params = {"alpha": -1, "beta": 1, "gamma": 0.3, "delta":0.15, "omega":1}
# and decomment the follow

# for k in range (1,100):
#     PoincareMap(sol1,t,params["omega"],1+3*k)



# ==============================================================================
#  PLOT
plt.figure(dpi=150)
plt.figure(figsize=(10,6))
title = 'Duffing Oscillator \n ' + r'$\ddot{x} +\delta \dot{x}+\alpha x +\beta x^3 = \gamma \cos(\omega t)$'
plt.title(title)
plt.plot(t,x_RK4, "r",label="RK4")
plt.plot(t,x_AB3, "b", label="AB3")
plt.xlabel("Time")
plt.ylabel("x")
plt.grid()
plt.legend()
plt.show()

import numpy as np
import matplotlib.pyplot as plt

# Define equation
def LorenzAttractor(x, params):
    rho = params["rho"]
    sigma = params["sigma"]
    beta = params["beta"]

    xdot = np.array([sigma*(x[1]-x[0]), rho*x[0]-x[1]-x[0]*x[2],-beta*x[2]+x[0]*x[1]])

    return xdot

def RungeKutta4(f, x0, t0, tf, dt):

    t = np.arange(t0, tf, dt)
    nt = t.size
    nx = x0.size

    x=np.zeros((nx,nt))
    print(x)

    x[:,0] = x0

    for k in range(nt -1):
        k1 = dt*f(t[k], x[:,k])
        k2 = dt*f(t[k] + dt/2, x[:,k]+ k1/2)
        k3 = dt*f(t[k] + dt/2, x[:,k]+ k2/2)
        k4 = dt*f(t[k] + dt, x[:,k]+ k3)

        dx = 1/6*(k1+2*k2+2*k3+k4)

        x[:,k+1] = x[:,k] + dx

    return x,t


params = {"rho": 28, "sigma": 10, "beta": 8/3}

f = lambda t,x : LorenzAttractor(x, params)

x0=np.array([1,-1,1])
t0=0
tf=100
dt=0.01

#Solve the equation

x,t = RungeKutta4(f,x0,t0,tf,dt)

# Plot results

plt.subplot(1,2,1)
plt.plot(t,x[0,:], "r",label="x")
plt.plot(t,x[1,:], "b", label="y")
plt.plot(t,x[2,:], "g", label="z")
plt.xlabel("Time")
plt.grid()
plt.legend()

plt.subplot(1,2,2)
plt.plot(x[0,:],x[1,:])
plt.xlabel("Preys")
plt.ylabel("Predators")
plt.grid()

plt.show()
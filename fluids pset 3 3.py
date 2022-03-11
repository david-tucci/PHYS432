"""
@author: David Tucci
@collab: Jules Faucher, Mattias Lazda 
2022-03-10
This code was based off prof.Lee notes
"""


import numpy as np
import matplotlib.pyplot as pl

# Parameters
time_step = 9
space_step = 1
Ngrid = 120
Nsteps = 30000
g = 0.001
angle = 10 # in the degrees
D = 0.1

beta = D*time_step/space_step**2

# Acceleration due to gravity
a =g*np.sin(angle*np.pi/180)

# Define variables 
x = np.arange(0, Ngrid*1., space_step) # multiplying by 1. to make sure this is an array of floats not integers
f = np.zeros(len(x))  

# Expected velocity determined in class lecture 11.8
u_analytic = -a/D *(0.5*np.power(x,2)-x[-1]*x) 

pl.ion()
fig, ax = pl.subplots(1,1)   


# Updating these plotting objects
plt, = ax.plot(x, f, 'ro')
ax.plot(x,u_analytic)

# # this draws the objects on the plot
fig.canvas.draw()

# Setting up matrices for diffusion operator
A1 = np.eye(Ngrid) * (1.0 + 2.0 * beta) + np.eye(Ngrid, k=1) * -beta + np.eye(Ngrid, k=-1) * -beta

# no slip boundary condition 
A1[0][0] = 1
A1[0][1] = 0
    
# no stress boundary condition
A1[Ngrid-1][Ngrid-1] = 1 + beta

ax.set_xlabel('distance')
ax.set_ylabel('velocity')

# Looping function f
for ct in range(Nsteps):
    # Solving for the next timestep
    f = np.linalg.solve(A1, f)
    
    # forces
    f[1:] += a*time_step
    
    # To make it quicker
    if ct%1000==0:
    # update the plot
        plt.set_ydata(f)
    
        fig.canvas.draw()
        pl.pause(0.001)
        ax.set_title('step = {0}'.format(ct))

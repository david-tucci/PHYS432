# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 17:35:54 2022

@author: David Tucci
@collab: Mattias Lazda, Jules Faucher 
"""

# Questions are answered at the bottom of this code and snapshots are uploaded on crowdmark


import numpy as np
import matplotlib.pyplot as pl

# Set up the grid, time and grid spacing, and the sound speed squared
Ngrid = 100
Nsteps = 5000
dt = 0.012
dx = 1.5
gamma = 5/3

x = np.arange(Ngrid) * dx # grid
f1 = np.ones(Ngrid) # rho
f2 = np.zeros(Ngrid) # rho x u
f3 = np.ones(Ngrid) # rho x etot
u = np.zeros(Ngrid+1) # advective velocity (keep the 1st and last element zero)

#Pressure
P = (gamma -1)/gamma * (f3 -0.5*f2**2/f1)

# Sound speed
cs2 = gamma * P/f1

# velocity and mach number
u_wave = f2/f1
mach = u_wave/np.sqrt(cs2)

def advection(f, u, dt, dx):
    # calculating flux terms
    J = np.zeros(len(f)+1) # keeping the first and the last term zero
    J[1:-1] = np.where(u[1:-1] > 0, f[:-1] * u[1:-1], f[1:] * u[1:-1])
    f = f - (dt / dx) * (J[1:] - J[:-1]) #update

    return f

# Apply initial Gaussian perturbation
Amp, sigma = 1000, Ngrid/10
f3 = f3 + Amp * np.exp(-(x - x.max()/2) ** 2 / sigma ** 2)

# plotting
pl.ion()
fig, ax = pl.subplots(2,1)

# Strong shock limit
rho2 = np.empty(Ngrid)
rho2.fill((gamma+1)/(gamma-1))


# Top graph density
x1, = ax[0].plot(x, f1, 'r-')
x11, = ax[0].plot(x, rho2, label = 'Strong Shock')
# ax[0].set_xlim([0, dx*Ngrid+1])
ax[0].set_ylim([0, 5])
ax[0].set_xlabel('x')
ax[0].set_ylabel('Density')
ax[0].legend()

# Bottom Mach number
x2, = ax[1].plot(x, mach, 'r-')
# ax[1].set_xlim([0, dx*Ngrid+1])
ax[1].set_ylim([-2, 2])
ax[1].set_xlabel('x')
ax[1].set_ylabel('Mach number')


fig.canvas.draw()

for ct in range(Nsteps):
    # advection velocity at the cell interface
    u[1:-1] = 0.5 * ((f2[:-1] / f1[:-1]) + (f2[1:] / f1[1:]))
    u[0] = 0
    u[-1] = 0
    
    # update density and momentum
    f1 = advection(f1, u, dt, dx)
    f2 = advection(f2, u, dt, dx)
    
    # Pressure
    P = (gamma-1)/gamma * (f3 -0.5*f2**2/f1)
    
    # Compute sound speed
    cs2 = gamma*P/f1
    
    # Pressure gradient
    f2[1:-1] = f2[1:-1] - 0.5*(dt/dx)* (P[2:] - P[:-2])
    
    # correct for source term at the boundary (reflective)
    f2[0] = f2[0] - 0.5 * (dt / dx) * (P[1] - P[0])
    f2[-1] = f2[-1] - 0.5 * (dt / dx) * (P[-1] - P[-2])

    # re-calculate advection velocities
    u[1:-1] = 0.5 * ((f2[:-1]/ f1[:-1]) + (f2[1:]/f1[1:]))
    
    # Advect energy
    f3 = advection(f3, u, dt, dx)
    
    # Recalculate pressure
    P = (gamma-1)/gamma * (f3 -0.5*f2**2/f1)
    
    # Define so easier to write
    u_wave = f2/f1
    Pu = P*u_wave
    
    # add the source term to energy
    f3[1:-1] = f3[1:-1] - 0.5*(dt/dx) * (Pu[2:] - Pu[:-2])
    
    # correct for source term at the boundary (reflective)
    f3[0] = f3[0] - 0.5 * (dt / dx) *  (Pu[1] - Pu[0])
    f3[-1] = f3[-1] - 0.5 * (dt / dx) * (Pu[-1] - Pu[-2])
    
    # Update Pressure, sound wave and mach number
    P = (gamma-1)/gamma * (f3 -0.5*f2**2/f1)
    cs2 = gamma*P/f1
    u_wave = f2/f1
    mach = u_wave/np.sqrt(cs2)    
    
    # update the plot
    x1.set_ydata(f1)
    x2.set_ydata(mach)
    x11.set_ydata(rho2)
    fig.canvas.draw()
    pl.pause(0.001)



## Question 1: Ratio of densities pre and post shock
# The ratio of densities is x = diff_rho = (gamma -1)/(gamma+1) + 2/(mach**2(1+gamma))
# ,but since we are in the strong shcok solution the second term vanishes and the
# ratio eaquals 1/4 for a monatomic gas. It does agree with what we got in class
# since the densities never surpass this limit, until the boundary which is to be
# expected since it "reflects" off. To see that it does not pass the limit and
# that it will at the boundaries two snapshots are attached.


# Question 2: Width of shock 
# The width of the shock is set by the viscosity divided by the velocity, since
# the mean free path divided by the mach number gave the width of the shock. The
# viscosity in numerical form is visc = dx^2/2*dt. Basically, the collisonal
# mean free path scales with viscosity and that's why the width of the shock
# also depends on the viscosity of the fluid. For dt = 0.012 and dx = 1.5, the
#  width of the shock is around 20 units as can be seen in the snapshot attached.
   
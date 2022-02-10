"""
@author: David Tucci
@collab: Jules Faucher, Mattias Lazda 
2022-02-10
"""
import numpy as np
import matplotlib.pyplot as pl

dt = 1  # Time step
Nsteps = 100

## Setting up initial conditions (vortex centres and circulations)
# Vortex rings

y_v = np.array([5,5,-5,-5],dtype = float) # y position of 4 vortices
x_v = np.array([-5,5,-5,5], dtype = float)
k_v = np.array([1,1,-1,-1],dtype = float)

# Setting up the plot
pl.ion()
fig, ax =pl.subplots(1,1)
# mark the initial positions of vortices

p, = ax.plot(x_v,y_v, 'k+', markersize = 10)

# draw initial velocity streamline
ngrid = 30 # insert the dimension of your simulation grid
Y, X = np.mgrid[-ngrid:ngrid:360j, -ngrid:ngrid:360j]
#360j sets the resolution of the cartesian grid ; play around with it
vel_x = np.zeros(np.shape(Y)) #this holds y-velocity
vel_y = np.zeros(np.shape(Y)) #this holds y-velocity

# masking radius for better visualization of the vortex centres
r_mask = 1 # insert the radius of the mask around the vortex centres
# within this mask, you will not plot any streamline
# so that you can see more clearly the movement of the vortex centres

for i in range(len(x_v)):#looping over each vortex
    # insert lines for computing total velocity field
    # Difference between X of the grid and position of vortex
    x = X-x_v[i] 
    y = Y-y_v[i]
    
    # Since velocity is for a line vortex u = k/r in the azimuthal direction
    r = np.sqrt(np.power(x,2)+np.power(y,2))
    
    xvelocity = (k_v[i]/r**2)*-y                                                              
    yvelocity = (k_v[i]/r**2)*x
    
    vel_x += xvelocity
    vel_y +=  yvelocity
     
    #insert lines for masking (set the masking area to NaN)
    vel_x[r <= r_mask] = np.nan
    vel_y[r <= r_mask] = np.nan

    
#set up the boundaries of the simulation box

ax.set_xlim([-ngrid, ngrid])
ax.set_ylim([-ngrid, ngrid])

#initial plot of the streamlines
ax.streamplot(X,Y, vel_x, vel_y, color='b', density = [1,1])
#play with density

fig.canvas.draw()



# Evolution
count = 0
while count < Nsteps:
    
    ## compute and update advection velocity
     
    
    # insert lines to re-initialize the total velocity field
    vx = np.zeros(len(x_v)) 
    vy = np.zeros(len(x_v))
    
    for i in range (len(x_v)):
        
         
        #insert lines to compute the total advection velocity on each vortex
        cur_x = x_v[i] # Current positions
        cur_y = y_v[i]
        filtre = np.ones(len(x_v),dtype = bool)
        filtre[i] = False
        # Other positions and k values of the other vortices 
        other_x = x_v[filtre] 
        other_y = y_v[filtre]
        other_k = k_v[filtre]
        
        
        x = cur_x - other_x
        y = cur_y - other_y
        
        r = np.sqrt(x*x+y*y)
        
        # Total velocity from other vortices
        vx[i] = np.sum((-other_k*y/r**2))                                                             
        vy[i] = np.sum((other_k*x/r**2))
    
    # insert lines to update the positions of vortices
    x_v += vx*dt
    y_v += vy*dt  
          
    #insert lines for masking (set the masking area to NaN)
    mask_filter = np.ones(X.shape,bool)
         
         
    # insert lines to re - initialize the total velocity field
    vel_x = np.zeros(np.shape(Y)) 
    vel_y = np.zeros(np.shape(Y))    


    for i in range (len(x_v)):
        
        #insert lines for computing total velocity field    
        x = X-x_v[i]
        y = Y-y_v[i]
        
        r = np.sqrt(np.power(x,2)+np.power(y,2))
        
        xvelocity = (k_v[i]/r**2)*-y                                                              
        yvelocity = (k_v[i]/r**2)*x
        
        vel_x += xvelocity
        vel_y +=  yvelocity
         
        #insert lines for masking (set the masking area to NaN)
        vel_x[r <= r_mask] = np.nan
        vel_y[r <= r_mask] = np.nan

    ## update plot
    # this clears out the previous streamlines
    ax.collections = []
    ax.patches = []
    
    p.set_xdata(x_v)
    p.set_ydata(y_v)
    
    ax.streamplot(X,Y, vel_x, vel_y, color='b', density = [1,1])
    
    fig.canvas.draw()
    
    pl.pause(0.001) # play with
    count +=1
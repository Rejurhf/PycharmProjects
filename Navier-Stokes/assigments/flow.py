from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np
import math as mth

class Flow1:

    def __init__(self):
        print("flow")


    def runSimulation(self):
        print("run")

        nx = 201  # x-points
        ny = 101  # y-points
        nit= 50
        c =  1                       # phase propagation speed
        x_span = 20.0
        y_span = 10.0
        dx = x_span/(nx-1)           # size of x-grid
        dy = y_span/(ny-1)           # size of y-grid
        x = np.linspace(0,x_span,nx) # last point included, so exactly nx points
        y = np.linspace(0,y_span,ny) # last point included, so exactly ny points
        X,Y = np.meshgrid(x,y)       # makes 2-dimensional mesh grid

        botb  = 40           # bottom boundary of obstacle
        dpth  = 20           # obstacle depth

        lftb   = 70          # left boundary of obstacle
        wdth   = 5           # obstacle width

        Re = 50              # range from 10s to 100s
        nt = 1000            # timesteps

        u_start = 1          # initial velocity at the start
        rho = 1              # density
        nu = ((dy*dpth)*u_start)/Re  # viscosity (UL/Re, Re = UL/nu, original value: 0.1)
        dt = 0.001               # timesteps

        qres = 4                 # quiver plot resolution

        v = np.zeros((ny, nx))
        u = np.ones((ny, nx))    # for u-velocity I initialise to 1 everywhere

        p = np.zeros((ny, nx))
        # b = np.zeros((ny, nx))

        # u, v, p = cavityFlow(nt, u, v, dt, dx, dy, p, rho, nu, botb, dpth, lftb, wdth, X, Y, u_start)

        # Plot the last figure on screen
        fig = plt.figure(figsize=(100,50), dpi=25)
        plt.contourf(X, Y, p, alpha=0.5)  # alpha - background intensity
        plt.tick_params(axis='both', which='major', labelsize=40)
        cbar = plt.colorbar()
        cbar.ax.tick_params(labelsize=50)
        # plt.contour(X, Y, p)
        plt.quiver(X[::qres,::qres],Y[::qres,::qres],u[::qres,::qres],v[::qres,::qres]) ##plotting velocity
        # plt.broken_barh([(x[lftb+1],x[lftb+wdth-2]-x[lftb+1])], (y[botb+1],y[botb+dpth-2]-y[botb+1]), hold=None, facecolors='grey', alpha=0.8)
        plt.xlabel('X')
        plt.ylabel('Y');
        plt.title('time_step = ' + str(nt) + ' nu = ' + str(nu), fontsize=40)
        plt.show()
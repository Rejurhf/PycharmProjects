import numpy
from matplotlib import pyplot
import time
import sys


class LinearConvection1D:
    '''
    Linear Convection 1D
    run_simulation(nx, nt, dt, c)
    nx - number of grid points; nt - number of timesteps we want to calculate;
    dt - amount of time each timestep covers (delta t); c - wave speed
    '''
    def __init__(self):
        print("Init Linear Convection 1D")

    def run_simulation(self, _nx=41, _nt=25, _dt=.025, _c=1):
        nx = _nx  # number of grid points
        dx = 2 / (nx - 1)  # distance between any pair of adjacent grid points
        nt = _nt  # nt is the number of timesteps we want to calculate
        dt = _dt  # dt is the amount of time each timestep covers (delta t)
        c = _c  # assume wave speed of c = 1

        u = numpy.ones(nx)  # numpy function ones()
        u[int(.5 / dx):int(1 / dx + 1)] = 2  # setting u = 2 between 0.5 and 1 as per our I.C.s
        print(u)

        pyplot.plot(numpy.linspace(0, 2, nx), u)  # Show array u
        pyplot.show()

        un = numpy.ones(nx)  # initialize a temporary array

        for n in range(0, nt):  # loop for values of n from 0 to nt, so it will run nt times
            un = u.copy()  # copy the existing values of u into un
            for i in range(1, nx):  # you can try commenting this line and...
                u[i] = un[i] - c * dt / dx * (un[i] - un[i - 1])

            if n % 2 == 0:
                pyplot.plot(numpy.linspace(0, 2, nx), u)
                pyplot.show()

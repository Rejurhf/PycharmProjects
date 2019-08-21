import numpy
from matplotlib import pyplot


class NonlinearConvection1D:
    '''
    Nonlinear Convection 1D
    '''
    def __init__(self):
        print("Init Nonlinear Convection 1D")

    def run_simulation(self, _nx=41, _nt=20, _dt=.025, _c=1):
        nx = _nx  # number of grid points
        nt = _nt  # nt is the number of time steps we want to calculate
        dt = _dt  # dt is the amount of time each time step covers (delta t)
        dx = 2 / (nx - 1)  # distance between any pair of adjacent grid points

        u = numpy.ones(nx)  # initialize u with every value equal to 1.
        u[int(.5 / dx): int(1 / dx + 1)] = 2  # set u = 2 between 0.5 and 1 as per our I.C.s

        un = numpy.ones(nx)  # initialize placeholder array un, to hold the time-stepped solution

        for n in range(0, nt):  # loop for values of n from 0 to nt, so it will run nt times
            un = u.copy()  # copy the existing values of u into un
            for i in range(1, nx):  # you can try commenting this line and...
                u[i] = un[i] - c * dt / dx * (un[i] - un[i - 1])

            if n % 2 == 0:
                pyplot.plot(numpy.linspace(0, 2, nx), u)
                pyplot.show()
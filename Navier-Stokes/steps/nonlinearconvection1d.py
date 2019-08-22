import numpy as np
from matplotlib import pyplot as plt
import time
import sys


class NonlinearConvection1D:
    """
    Nonlinear Convection 1D
    """

    iteration_number = 0

    def __init__(self):
        print("Init Linear Convection 1D")

    def run_simulation(self, _grid_points=41, _grid_length=10, _nt=400, _dt=.025):
        grid_points = _grid_points  # number of grid points
        grid_length = _grid_length  # length of grid
        dx = grid_length / (grid_points - 1)  # distance between any pair of adjacent grid points
        nt = _nt  # nt is the number of time steps we want to calculate
        dt = _dt  # dt is the amount of time each time step covers (delta t)

        u = np.ones(grid_points)  # numpy function ones()
        u[int(.5 / dx):int(1 / dx + 1)] = 2  # setting u = 2 between 0.5 and 1 as per our I.C.s

        self.draw_1d_plot(grid_length, grid_points, u, '1D Linear Convection t=0')

        un = np.ones(grid_points)

        for n in range(nt):  # Runs however many timesteps you set earlier
            un = u.copy()  # copy the u array to not overwrite values
            for i in range(1, grid_points):
                u[i] = un[i] - un[i] * dt / dx * (un[i] - un[i - 1])

            if n % 20 == 0:
                self.draw_1d_plot(grid_length, grid_points, u, '1D Linear Convection t={}'.format(n/40))


    def draw_1d_plot(self, grid_length, grid_points, u, text):
        self.iteration_number += 1
        print("\n", self.iteration_number, ". --------------------------------------------")
        print("Sum of points in u: ", sum(u))
        print(u)
        plt.plot(np.linspace(0, grid_length, grid_points), u)  # Show array u
        plt.ylim(1, 2.1)
        plt.xlabel('x')
        plt.ylabel('u')
        plt.title(text)
        plt.show()

import numpy as np
from matplotlib import pyplot as plt
import sympy


class BurgersEquation1D:
    """
    Burgers' Equation 1D
    """

    iteration_number = 0

    def __init__(self):
        print("Init Burgers' Equation 1D")

    def run_simulation(self, _grid_points=101, _grid_length=2, _nt=150, _nu=.07):
        # Symbolic calculations
        x, nu, t = sympy.symbols('x nu t')
        phi = (sympy.exp(-(x - 4 * t) ** 2 / (4 * nu * (t + 1))) +
               sympy.exp(-(x - 4 * t - 2 * np.pi) ** 2 / (4 * nu * (t + 1))))
        phiprime = phi.diff(x)
        u = -2 * nu * (phiprime / phi) + 4
        ufunc = sympy.utilities.lambdify((t, x, nu), u)

        grid_points = _grid_points  # number of grid points
        grid_length = _grid_length  # length of grid
        nt = _nt  # nt is the number of time steps we want to calculate
        nu = _nu  # viscosity of the system
        dx = grid_length * np.pi / (grid_points - 1)  # distance between any pair of adjacent grid points
        dt = dx * nu  # Dynamically scaling dt based on grid size to ensure convergence
        print("dt = ", dt)

        # Initializing the array containing the shape of our initial conditions
        x = np.linspace(0, 2 * np.pi, grid_points)
        un = np.empty(grid_points)
        t = 0

        u = np.asarray([ufunc(t, x0, nu) for x0 in x])
        u_anal = np.asarray([ufunc(t * dt, xi, nu) for xi in x])
        self.draw_1d_plot(x, u, u_anal, "1D Burgers' Equation t=0")

        for n in range(nt):  # Runs however many timesteps you set earlier
            un = u.copy()  # copy the u array to not overwrite values
            for i in range(1, grid_points - 1):
                u[i] = un[i] - un[i] * dt / dx * (un[i] - un[i - 1]) + nu * (dt / dx ** 2) * (
                            un[i + 1] - 2 * un[i] + un[i - 1])

            u[0] = un[0] - un[0] * dt / dx * (un[0] - un[-2]) + nu * (dt / dx ** 2) * (un[1] - 2 * un[0] + un[-2])
            u[-1] = u[0]

            if n%15 == 0:
                u_anal = np.asarray([ufunc(n * dt, xi, nu) for xi in x])
                self.draw_1d_plot(x, u, u_anal, "1D Burgers' Equation t={}".format(n/15 + 1))





    def draw_1d_plot(self, x, u, u_anal, text):
        self.iteration_number += 1
        print("\n{}. --------------------------------------------".format(self.iteration_number))
        print("Sum of points in u: ", sum(u))
        print(u)

        plt.plot(x, u, marker='o', lw=2, label='Computational')
        plt.plot(x, u_anal, label='Analytical')
        plt.ylim(0, 8)
        plt.xlabel('x')
        plt.ylabel('u')
        plt.title(text)
        plt.show()
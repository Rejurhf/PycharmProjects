import numpy as np
from matplotlib import pyplot
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D


class NonlinearConvection2D:
    """
    Nonlinear Convection 2D
    """

    iteration_number = 0

    def __init__(self):
        print("Init Nonlinear Convection 2D")


    def run_simulation(self, _nx=21, _ny=21, _grid_len=2, _nt=50, _c=1, _sigma=.2):
        nx = _nx
        ny = _ny
        grid_len = _grid_len
        nt = _nt  # nt is the number of time steps we want to calculate
        c = _c  # assume wave speed of c = 1
        sigma = _sigma
        dx = grid_len / (nx - 1)  # distance between any pair of adjacent grid points
        dy = grid_len / (ny - 1)  # distance between any pair of adjacent grid points
        dt = sigma * dx  # dt is the amount of time each time step covers (delta t)

        x = np.linspace(0, grid_len, nx)
        y = np.linspace(0, grid_len, ny)

        u = np.ones((ny, nx))  #create a 1xn vector of 1's
        v = np.ones((ny, nx))

        # Initiallizing the array containing the shape of our initial conditions
        u[int(.5 / dy):int(1 / dy + 1), int(.5 / dx):int(1 / dx + 1)] = 2
        v[int(.5 / dy):int(1 / dy + 1), int(.5 / dx):int(1 / dx + 1)] = 2

        self.draw_3d_plot(x, y, u, "2D Linear Convection at t=0")

        for n in range(nt + 1):  ##loop across number of time steps
            un = u.copy()
            vn = u.copy()
            row, col = u.shape
            for j in range(1, row):
                for i in range(1, col):
                    u[1:, 1:] = (un[1:, 1:] -
                                 (un[1:, 1:] * c * dt / dx * (un[1:, 1:] - un[1:, :-1])) -
                                 vn[1:, 1:] * c * dt / dy * (un[1:, 1:] - un[:-1, 1:]))

                    v[1:, 1:] = (vn[1:, 1:] -
                                 (un[1:, 1:] * c * dt / dx * (vn[1:, 1:] - vn[1:, :-1])) -
                                 vn[1:, 1:] * c * dt / dy * (vn[1:, 1:] - vn[:-1, 1:]))

                    u[0, :] = 1
                    u[-1, :] = 1
                    u[:, 0] = 1
                    u[:, -1] = 1

                    v[0, :] = 1
                    v[-1, :] = 1
                    v[:, 0] = 1
                    v[:, -1] = 1
            if n%10 == 0 and n > 0:
                self.draw_3d_plot(x, y, u, "2D Linear Convection at t={}".format(n/10))


    def draw_3d_plot(self, x, y, u, text):
        self.iteration_number += 1
        print("\n", self.iteration_number, ". ---------------------------------------------")
        print(u)

        fig = pyplot.figure(figsize=(11, 7), dpi=100)
        ax = fig.gca(projection='3d')
        X, Y = np.meshgrid(x, y)
        surf = ax.plot_surface(X, Y, u[:], cmap=cm.viridis)
        ax.set_xlabel('$x$')
        ax.set_ylabel('$y$')
        ax.set_zlabel('$u$')
        ax.set_zlim(1, 2)
        ax.text2D(0.35, 0.95, text, transform=ax.transAxes)
        pyplot.show()



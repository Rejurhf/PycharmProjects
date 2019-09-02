import numpy
from matplotlib import pyplot, cm
from mpl_toolkits.mplot3d import Axes3D


class BurgersEquation2D:
    """
    Burgers' Equation 2D
    """

    iteration_number = 0

    def __init__(self):
        print("Init Burgers' Equation 2D")

    def run_simulation(self, _nx=41, _ny=41, _grid_len=2, _nt=3840, _nu=.01, _c=1, _sigma=.0009):
        nx = _nx
        ny = _ny
        grid_len = _grid_len
        nt = _nt  # nt is the number of time steps we want to calculate
        nu = _nu
        c = _c  # assume wave speed of c = 1
        sigma = _sigma
        dx = grid_len / (nx - 1)  # distance between any pair of adjacent grid points
        dy = grid_len / (ny - 1)  # distance between any pair of adjacent grid points
        dt = sigma * dx * dy / nu  # dt is the amount of time each time step covers (delta t)

        x = numpy.linspace(0, grid_len, nx)
        y = numpy.linspace(0, grid_len, ny)

        u = numpy.ones((ny, nx))  # create a 1xn vector of 1's
        v = numpy.ones((ny, nx))
        comb = numpy.ones((ny, nx))

        # Assign initial conditions
        # set hat function I.C. : u(.5<=x<=1 && .5<=y<=1 ) is 2
        u[int(.5 / dy):int(1 / dy + 1), int(.5 / dx):int(1 / dx + 1)] = 2
        v[int(.5 / dy):int(1 / dy + 1), int(.5 / dx):int(1 / dx + 1)] = 2

        self.draw_3d_plot(x, y, u, v, "2D Burgers' Equation t=0")

        for n in range(nt + 1):  ##loop across number of time steps
            un = u.copy()
            vn = v.copy()

            u[1:-1, 1:-1] = (un[1:-1, 1:-1] -
                             dt / dx * un[1:-1, 1:-1] *
                             (un[1:-1, 1:-1] - un[1:-1, 0:-2]) -
                             dt / dy * vn[1:-1, 1:-1] *
                             (un[1:-1, 1:-1] - un[0:-2, 1:-1]) +
                             nu * dt / dx ** 2 *
                             (un[1:-1, 2:] - 2 * un[1:-1, 1:-1] + un[1:-1, 0:-2]) +
                             nu * dt / dy ** 2 *
                             (un[2:, 1:-1] - 2 * un[1:-1, 1:-1] + un[0:-2, 1:-1]))

            v[1:-1, 1:-1] = (vn[1:-1, 1:-1] -
                             dt / dx * un[1:-1, 1:-1] *
                             (vn[1:-1, 1:-1] - vn[1:-1, 0:-2]) -
                             dt / dy * vn[1:-1, 1:-1] *
                             (vn[1:-1, 1:-1] - vn[0:-2, 1:-1]) +
                             nu * dt / dx ** 2 *
                             (vn[1:-1, 2:] - 2 * vn[1:-1, 1:-1] + vn[1:-1, 0:-2]) +
                             nu * dt / dy ** 2 *
                             (vn[2:, 1:-1] - 2 * vn[1:-1, 1:-1] + vn[0:-2, 1:-1]))

            u[0, :] = 1
            u[-1, :] = 1
            u[:, 0] = 1
            u[:, -1] = 1

            v[0, :] = 1
            v[-1, :] = 1
            v[:, 0] = 1
            v[:, -1] = 1

            if n%384 == 0 and n > 0:
                self.draw_3d_plot(x, y, u, v, "2D Burgers' Equation t={}".format(n/12))


    def draw_3d_plot(self, x, y, u, v, text):
        self.iteration_number += 1
        print("\n", self.iteration_number, ". ---------------------------------------------")
        print(u)

        fig = pyplot.figure(figsize=(11, 7), dpi=100)
        ax = fig.gca(projection='3d')
        X, Y = numpy.meshgrid(x, y)
        ax.plot_surface(X, Y, u[:], cmap=cm.viridis, rstride=1, cstride=1)
        ax.plot_surface(X, Y, v[:], cmap=cm.viridis, rstride=1, cstride=1)
        ax.set_xlabel('$x$')
        ax.set_ylabel('$y$')
        ax.set_zlabel('$u$')
        ax.set_zlim(1, 2)
        ax.text2D(0.35, 0.95, text, transform=ax.transAxes)
        pyplot.show()
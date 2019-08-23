import numpy
from matplotlib import pyplot, cm
from mpl_toolkits.mplot3d import Axes3D


class LaplaceEquation2D:
    """
    Laplace Equation 2D
    """

    iteration_number = 0

    def __init__(self):
        print("Init Laplace Equation 2D")

    def run_simulation(self, _nx=31, _ny=31, _nt=200, _grid_len=2, _c=1):
        nx = _nx
        ny = _ny
        nt = _nt
        grid_len = _grid_len
        c = _c  # assume wave speed of c = 1
        dx = grid_len / (nx - 1)  # distance between any pair of adjacent grid points
        dy = grid_len / (ny - 1)  # distance between any pair of adjacent grid points

        ##initial conditions
        p = numpy.zeros((ny, nx))  # create a XxY vector of 0's

        ##plotting aids
        x = numpy.linspace(0, 2, nx)
        y = numpy.linspace(0, 1, ny)

        ##boundary conditions
        p[:, 0] = 0  # p = 0 @ x = 0
        p[:, -1] = y  # p = y @ x = 2
        p[0, :] = p[1, :]  # dp/dy = 0 @ y = 0
        p[-1, :] = p[-2, :]  # dp/dy = 0 @ y = 1

        self.plot2d(x, y, p, "2D Laplace Equation t=0")

        for n in range(nt):
            p = self.laplace2d(p, y, dx, dy)

            if n%20 == 0 and n > 0:
                self.plot2d(x, y, p, "2D Laplace Equation t={}".format(n/10))


    def laplace2d(self, p, y, dx, dy):
        pn = numpy.empty_like(p)
        pn = p.copy()
        p[1:-1, 1:-1] = ((dy ** 2 * (pn[1:-1, 2:] + pn[1:-1, 0:-2]) +
                          dx ** 2 * (pn[2:, 1:-1] + pn[0:-2, 1:-1])) /
                         (2 * (dx ** 2 + dy ** 2)))

        p[:, 0] = 0  # p = 0 @ x = 0
        p[:, -1] = y  # p = y @ x = 2
        p[0, :] = p[1, :]  # dp/dy = 0 @ y = 0
        p[-1, :] = p[-2, :]  # dp/dy = 0 @ y = 1

        return p


    def plot2d(self, x, y, p, text):
        self.iteration_number += 1
        print("\n", self.iteration_number, ". ---------------------------------------------")
        print(p)

        fig = pyplot.figure(figsize=(11, 7), dpi=100)
        ax = fig.gca(projection='3d')
        X, Y = numpy.meshgrid(x, y)
        surf = ax.plot_surface(X, Y, p[:], rstride=1, cstride=1, cmap=cm.viridis,
            linewidth=0, antialiased=False)
        ax.set_xlabel('$x$')
        ax.set_ylabel('$y$')
        ax.set_zlabel('$u$')
        ax.set_xlim(0, 2)
        ax.set_ylim(0, 1)
        ax.set_zlim(0, 1)
        ax.view_init(30, 225) # turn plot
        ax.text2D(0.35, 0.95, text, transform=ax.transAxes)
        pyplot.show()
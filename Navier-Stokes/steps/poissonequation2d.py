import numpy
from matplotlib import pyplot, cm
from mpl_toolkits.mplot3d import Axes3D


class PoissonEquation2D:
    """
    Poisson Equation 2D
    """

    iteration_number = 0

    def __init__(self):
        print("Init Poisson Equation 2D")

    def run_simulation(self, _nx=50, _ny=50, _nt=100, _xmin=0, _xmax=2, _ymin=0, _ymax=1):
        nx = _nx
        ny = _ny
        nt = _nt  # nt is the number of time steps we want to calculate
        xmin = _xmin
        xmax = _xmax
        ymin = _ymin
        ymax = _ymax

        dx = (xmax - xmin) / (nx - 1)  # distance between any pair of adjacent grid points
        dy = (ymax - ymin) / (ny - 1)  # distance between any pair of adjacent grid points

        # Initialization
        p = numpy.zeros((ny, nx))
        pd = numpy.zeros((ny, nx))
        b = numpy.zeros((ny, nx))
        x = numpy.linspace(xmin, xmax, nx)
        y = numpy.linspace(xmin, xmax, ny)

        # Source
        b[int(ny / 4), int(nx / 4)] = 100
        b[int(3 * ny / 4), int(3 * nx / 4)] = -100

        self.plot2d(x, y, p, "2D Laplace Equation t=0")

        for it in range(nt):
            pd = p.copy()

            p[1:-1, 1:-1] = (((pd[1:-1, 2:] + pd[1:-1, :-2]) * dy ** 2 +
                              (pd[2:, 1:-1] + pd[:-2, 1:-1]) * dx ** 2 -
                              b[1:-1, 1:-1] * dx ** 2 * dy ** 2) /
                             (2 * (dx ** 2 + dy ** 2)))

            p[0, :] = 0
            p[ny - 1, :] = 0
            p[:, 0] = 0
            p[:, nx - 1] = 0

            if it%10 == 0 and it > 0:
                self.plot2d(x, y, p, "2D Laplace Equation t={}".format(it/10))


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
        ax.set_zlim(-0.04, 0.04)
        ax.view_init(30, 225) # turn plot
        ax.text2D(0.35, 0.95, text, transform=ax.transAxes)
        pyplot.show()
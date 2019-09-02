import numpy
from matplotlib import pyplot, cm
from mpl_toolkits.mplot3d import Axes3D


class CavityFlow:
    """
    Cavity Flow
    """

    iteration_number = 0

    def __init__(self):
        print("Init Cavity Flow")

    def run_simulation(self, _nx=41, _ny=41, _nt=1000, _nit=50, _c=5):
        nx = _nx
        ny = _ny
        nt = _nt  # nt is the number of time steps we want to calculate
        nit = _nit
        c = _c

        dx = 2 / (nx - 1)
        dy = 2 / (ny - 1)
        x = numpy.linspace(0, 2, nx)
        y = numpy.linspace(0, 2, ny)
        X, Y = numpy.meshgrid(x, y)

        rho = 1.5
        nu = .09
        dt = .001

        u = numpy.zeros((ny, nx))
        v = numpy.zeros((ny, nx))
        p = numpy.zeros((ny, nx))
        b = numpy.zeros((ny, nx))

        # self.plot2d(X, Y, u, v, p, "Cavity Flow t=0")

        self.cavity_flow(nt, u, v, dt, dx, dy, p, rho, nu, nx, ny, nit, X, Y)


    def cavity_flow(self, nt, u, v, dt, dx, dy, p, rho, nu, nx, ny, nit, X, Y):
        un = numpy.empty_like(u)
        vn = numpy.empty_like(v)
        b = numpy.zeros((ny, nx))

        for n in range(nt):
            un = u.copy()
            vn = v.copy()

            b = self.build_up_b(b, rho, dt, u, v, dx, dy)
            p = self.pressure_poisson(p, dx, dy, b, nit)

            u[1:-1, 1:-1] = (un[1:-1, 1:-1] -
                             un[1:-1, 1:-1] * dt / dx *
                             (un[1:-1, 1:-1] - un[1:-1, 0:-2]) -
                             vn[1:-1, 1:-1] * dt / dy *
                             (un[1:-1, 1:-1] - un[0:-2, 1:-1]) -
                             dt / (2 * rho * dx) * (p[1:-1, 2:] - p[1:-1, 0:-2]) +
                             nu * (dt / dx ** 2 *
                                   (un[1:-1, 2:] - 2 * un[1:-1, 1:-1] + un[1:-1, 0:-2]) +
                                   dt / dy ** 2 *
                                   (un[2:, 1:-1] - 2 * un[1:-1, 1:-1] + un[0:-2, 1:-1])))

            v[1:-1, 1:-1] = (vn[1:-1, 1:-1] -
                             un[1:-1, 1:-1] * dt / dx *
                             (vn[1:-1, 1:-1] - vn[1:-1, 0:-2]) -
                             vn[1:-1, 1:-1] * dt / dy *
                             (vn[1:-1, 1:-1] - vn[0:-2, 1:-1]) -
                             dt / (2 * rho * dy) * (p[2:, 1:-1] - p[0:-2, 1:-1]) +
                             nu * (dt / dx ** 2 *
                                   (vn[1:-1, 2:] - 2 * vn[1:-1, 1:-1] + vn[1:-1, 0:-2]) +
                                   dt / dy ** 2 *
                                   (vn[2:, 1:-1] - 2 * vn[1:-1, 1:-1] + vn[0:-2, 1:-1])))

            u[0, :] = 0
            u[:, 0] = 0
            u[:, -1] = 0
            u[-1, :] = 1  # set velocity on cavity lid equal to 1
            v[0, :] = 0
            v[-1, :] = 0
            v[:, 0] = 0
            v[:, -1] = 0

            if n == 2:
                self.plot2d(X, Y, u, v, p, "Cavity Flow t=0")
            if n%100 == 0 and n > 0:
                self.plot2d(X, Y, u, v, p, "Cavity Flow t={}".format(n/20))




    def build_up_b(serlf, b, rho, dt, u, v, dx, dy):
        b[1:-1, 1:-1] = (rho * (1 / dt *
                                ((u[1:-1, 2:] - u[1:-1, 0:-2]) /
                                 (2 * dx) + (v[2:, 1:-1] - v[0:-2, 1:-1]) / (2 * dy)) -
                                ((u[1:-1, 2:] - u[1:-1, 0:-2]) / (2 * dx)) ** 2 -
                                2 * ((u[2:, 1:-1] - u[0:-2, 1:-1]) / (2 * dy) *
                                     (v[1:-1, 2:] - v[1:-1, 0:-2]) / (2 * dx)) -
                                ((v[2:, 1:-1] - v[0:-2, 1:-1]) / (2 * dy)) ** 2))

        return b


    def pressure_poisson(self, p, dx, dy, b, nit):
        pn = numpy.empty_like(p)
        pn = p.copy()

        for q in range(nit):
            pn = p.copy()
            p[1:-1, 1:-1] = (((pn[1:-1, 2:] + pn[1:-1, 0:-2]) * dy ** 2 +
                              (pn[2:, 1:-1] + pn[0:-2, 1:-1]) * dx ** 2) /
                             (2 * (dx ** 2 + dy ** 2)) -
                             dx ** 2 * dy ** 2 / (2 * (dx ** 2 + dy ** 2)) *
                             b[1:-1, 1:-1])

            p[:, -1] = p[:, -2]  # dp/dx = 0 at x = 2
            p[0, :] = p[1, :]  # dp/dy = 0 at y = 0
            p[:, 0] = p[:, 1]  # dp/dx = 0 at x = 0
            p[-1, :] = 0  # p = 0 at y = 2

        return p


    def plot2d(self, X, Y, u, v, p, text):
        self.iteration_number += 1
        print("\n", self.iteration_number, ". ---------------------------------------------")
        print(p)

        fig = pyplot.figure(figsize=(11, 7), dpi=100)
        # plotting the pressure field as a contour
        pyplot.contourf(X, Y, p, alpha=0.5, cmap=cm.viridis)
        pyplot.colorbar()
        # plotting the pressure field outlines
        pyplot.contour(X, Y, p, cmap=cm.viridis)
        # plotting velocity field
        pyplot.quiver(X[::2, ::2], Y[::2, ::2], u[::2, ::2], v[::2, ::2])
        pyplot.xlabel('X')
        pyplot.ylabel('Y')
        pyplot.title(text)
        pyplot.show()


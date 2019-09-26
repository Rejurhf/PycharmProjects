"""
# Use proper import for project wanted
# 1. 1D Linear Convection
from steps.linearconvection1d import LinearConvection1D
from steps.nonlinearconvection1d import NonlinearConvection1D
from steps.diffusion1d import Diffusion1D
from steps.burgersequation1d import BurgersEquation1D
from steps.linearconvection2d import LinearConvection2D
from steps.nonlinearconvection2d import NonlinearConvection2D
from steps.diffusion2d import Diffusion2D
from steps.burgersequation2d import BurgersEquation2D
from steps.laplaceequation2d import LaplaceEquation2D
from steps.poissonequation2d import PoissonEquation2D
from steps.cavityflow import CavityFlow
from steps.chanelflow import ChanelFlow

# 1. Linear Convection 1D
'''
run_simulation(grid_points, grid_length, nt, dt, c)
grid_points - number of grid points; nt - number of time steps we want to calculate;
dt - amount of time each tim estep covers (delta t); c - wave speed
'''
# LinearConvection1D().run_simulation()

# 2. Nonlinear Convection 1D
'''
run_simulation(grid_points, grid_length, nt, dt)
grid_points - number of grid points; nt - number of time steps we want to calculate;
dt - amount of time each time step covers (delta t)
'''
# NonlinearConvection1D().run_simulation()

# 3. Diffusion Equation 1D
# Diffusion1D().run_simulation()

# 4. Burgers' Equation 1D
# BurgersEquation1D().run_simulation()

# 5. Linear Convection 2D
# LinearConvection2D().run_simulation()

# 6. Nonlinear Convection 2D
# NonlinearConvection2D().run_simulation()

# 7. Diffusion Equation 2D
# Diffusion2D().run_simulation()

# 8. Burgers' Equation 2D
# BurgersEquation2D().run_simulation()

# 9. Laplace Equation 2D
# LaplaceEquation2D().run_simulation()

# 10. Poisson Equation 2D
# PoissonEquation2D().run_simulation()

# 11. Poisson Equation 2D
# CavityFlow().run_simulation()

# 12. Poisson Equation 2D
# ChanelFlow().run_simulation()

"""

from assigments.flow import Flow1
from assigments.flow2d import Flow2d


# Flow1().runSimulation()
Flow2d().runSimulation()


# Animation().runAnimation()
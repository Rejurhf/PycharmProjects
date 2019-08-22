
# Use proper import for project wanted
# 1. 1D Linear Convection
from steps.linearconvection1d import LinearConvection1D
from steps.nonlinearconvection1d import NonlinearConvection1D
from steps.diffusion1d import Diffusion1D

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


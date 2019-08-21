
# Use proper import for project wanted
# 1. 1D Linear Convection
from steps.linearconvection1d import LinearConvection1D
from steps.nonlinearconvection1d import NonlinearConvection1D

# 1. Linear Convection 1D
'''
run_simulation(nx, nt, dt, c)
nx - number of grid points; nt - number of timesteps we want to calculate;
dt - amount of time each timestep covers (delta t); c - wave speed
'''
LinearConvection1D().run_simulation()

# 2. Nonlinear Convection 1D
# NonlinearConvection1D.run_simulation()


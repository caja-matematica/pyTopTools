#!/usr/bin/python
# Implements the stochastic Hopf-bifurcation model illustrated in
# Week 308 of John Baez's This Week's Finds
# http://johncarlosbaez.wordpress.com/2010/12/24/this-weeks-finds-week-308/

import sys, random             # Python ships with these
import scipy, pylab     # these are extra

def hopf(x, y, beta, dt, lamb):
    """
    Compute the change in coordinates given the current position,
    the parameters which govern the stochastic Hopf dynamics and the 
    Euler integration timestep.
    """
    rsquared = x * x + y * y
    dx = (-y + beta * x - x * rsquared) * dt
    dy = (x + beta * y - y * rsquared) * dt
    if lamb != 0.0:
        sigma = sqrt(dt)
        dx += lamb * gauss(mu = 0.0, sigma = sigma)
        dy += lamb * gauss(mu = 0.0, sigma = sigma)
    return dx, dy

def split_param(param_string):
    """
    Split an item taken from the command line into a variable name and
    a value.
    """
    split = param_string.split("=")
    try:
        param_name = split[0].lower()
        param_value = float(split[1])
    except:
        param_name = param_string.lower()
        param_value = 0.0
    return (param_name, param_value)

def parse_command_line(command_line):
    """
    Use the command line arguments to specify the values of the
    simulation parameters; assign reasonable default values to
    all parameters left unspecified.
    """
    pairs = dict([split_param(item)
                  for item in command_line])
    try:
        beta = pairs["beta"]        # x-y coupling coefficient
    except:
        beta = 0.1
    try:
        lamb = pairs["lambda"]      # noise strength
    except:
        lamb = 0.0
    try:
        dt = pairs["dt"]            # timestep
    except:
        dt = 0.001
    try:
        x_init = pairs["x"]         # initial value for x
    except:
        x_init = 0.0
    try:
        y_init = pairs["y"]         # initial value for y
    except:
        y_init = 0.0
    try:
        T = pairs["T"]              # duration of the simulation
    except:
        T = int(1e4)
    return beta, lamb, dt, T, x_init, y_init

# pull functions out of libraries for our later convenience
gauss = random.gauss
sqrt = scipy.sqrt

# get command-line parameters
beta, lamb, dt, T, x_init, y_init  = parse_command_line(sys.argv[1:])

# initialize the arrays of coordinate values
aX = scipy.zeros(T)
aY = scipy.zeros(T)
aX[0] = x_init
aY[0] = y_init

# MAIN LOOP
for time_step in xrange(1, T):
    x = aX[time_step - 1]
    y = aY[time_step -1]
    dx, dy = hopf(x, y, beta, dt, lamb)
    aX[time_step] = x + dx
    aY[time_step] = y + dy

# display output
# Fig. 1: X vs. Y
pylab.plot(aX, aY)
pylab.xlabel("X")
pylab.ylabel("Y")
pylab.title("beta = " + str(beta)
            + ", lambda = " + str(lamb))
# Fig. 2: X as a function of time
pylab.figure()
pylab.plot(aX)
pylab.xlabel("Time")
pylab.ylabel("X")
pylab.title("beta = " + str(beta)
            + ", lambda = " + str(lamb))
pylab.show()

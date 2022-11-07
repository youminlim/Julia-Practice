# 2 Body Problem Solver

# Add Packages and Libraries
cd(@__DIR__)
using Pkg
Pkg.activate(".")

using Plots
include("myFunctions.jl")

# System parameters
G = 6.67 * 10 ^-11
m1 = 10*10^26
m2 = 10*10^26
y0 = [0, 0, 0, 3000*10^3, 0, 0, 10000, 20000, 30000, 0, 40000, 0] # In the form y0 = [x01 y01 z01 x02 y02 z02 dot(x01 y01 z01 x02 y02 z02)]
r = eulerian(y0[1:3], y0[4:6])

rk = 4 # Choice between 1, 2, 3 and 4

# Time interval
t0 = 0      # Initial time
tf = 480    # Final time
n = 480

t, y = rkSolver(rk, twoBody3D, t0, tf, y0, n)

plot2Body(t, y)

accelChecker(G, m1, m2, y)
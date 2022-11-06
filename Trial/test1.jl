# Getting started with Julia for experienced programmers - Notes

# This allows the compiler to find the location of the file within the directory
cd(@__DIR__)
# Calls the package manager in Julia to be used
using Pkg
# Uses the package manager to activate the packages within the environment/project this file is inside
Pkg.activate(".")

# Calls the packages to the file
using ForwardDiff, FiniteDiff

# Test code
f(x) = 2x^2 + x
ForwardDiff.derivative(f,2.0)
FiniteDiff.finite_difference_derivative(f,2.0)

# Creates packages from the project
using PkgTemplates
t = Template(user = "youminlim")
t("CFD_coursework")
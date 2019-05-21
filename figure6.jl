
include("ks.jl")

using Plots
using PyPlot
pygui(true)
pyplot()

# evaluate the method for different timesteps
# here, we take the timesteps 1/2^x for all x such that it is greater than 1/1000
hh = 2.0
hs = [1.0/hh]
while hh < 1000
    global hh
    hh *= 2
    push!(hs, 1/hh)
end
# compute the solutions for the different timesteps
etdrk4seq = [ks_etdrk4(hhh) for hhh in hs]

# we use the smallest timestep to approximate the exact solutions
exact_solution = etdrk4seq[10]
# the solutions we want to check convergence for
test_solutions = etdrk4seq[1:9]

# the time we want to test at (101 means time 150)
test_time = 101

# compute the infinity-norm for all test solutions
errors = []
for candidate in test_solutions
    global errors
    maxval = 0
    maxpack = 0
    for i in 1:128
        maxval = max(abs(exact_solution[2][test_time][i] - candidate[2][test_time][i]),maxval)
        maxpack = max(maxpack, abs(exact_solution[2][test_time][i]))
    end
    push!(errors, maxval/maxpack)
end

Plots.plot(hs[1:9], errors, xaxis=:log, yaxis=:log, markershape=:hexagon)

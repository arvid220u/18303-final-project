
include("kdv.jl")

using Plots
using PyPlot
pyplot()
pygui(true)

h = 1/1000
l = 2
N = 1024
# solve it! (from kdv.jl)
tt, uu = kdv_etdrk4(h, 1)

function extract_point(x,t)
    global h, l, N
    return uu[Int64(floor(t / h + 1))][Int64(floor(x / l  * N))]
end

l = 2
space_interval = (0.0,l)
N = 1024
x = space_interval[1] .+ space_interval[2] * (1:N)/N
t_series = tt[1:10:size(tt)[1]]#[t * 1 for t in 1:50]#size(uu)[1]]
x_series = x[1:10:size(x)[1]]#[x * 1 for x in 1:128]#N]


display(Plots.plot(x_series, t_series, extract_point, st=:surface, xlabel = "Space", ylabel="Time"))

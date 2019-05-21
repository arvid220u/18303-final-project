# solution of Kuramoto-Sivashinsky

include("ks.jl")

using PyPlot
using Plots
pyplot()
pygui(true)

uu = ks_etdrk4(1/4)[2]
N = 128

function extract_point(x,y)
    return uu[x][y]
end

x_series = [x * 1 for x in 1:size(uu)[1]]
y_series = [y * 1 for y in 1:N]

display(Plots.plot(x_series, y_series, extract_point, st=:surface, xticks = nothing, yticks = nothing))

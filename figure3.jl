
include("kdv.jl")

etdrk4 = kdv_etdrk4(1/1000, 1)

package = kdv_package()

space_interval = (0.0,2)
N = 128
x = space_interval[1] .+ space_interval[2] * (1:N)/N
x_series = x_series = x[1:size(x)[1]]

function etdrk4_value(x)
    ind = Int64(round(x * 128 / 2))
    return etdrk4[2][1001][8:8:1024][ind]
end
function package_value(x)
    ind = Int64(round(x * 128 / 2))
    return package[2][11][ind]
end

using Plots
using PyPlot
pyplot()
pygui(true)
Plots.plot(x_series, [etdrk4_value, package_value])

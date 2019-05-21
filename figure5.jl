
include("kdv.jl")

etdrk4 = kdv_etdrk4(1/500, 1)
etdrk4_small_contour = kdv_etdrk4(1/500, 0.000001)

package = kdv_package()

space_interval = (0.0,2)
N = 128
x = space_interval[1] .+ space_interval[2] * (1:N)/N
x_series = x_series = x[1:size(x)[1]]

function etdrk4_value(x)
    ind = Int64(round(x * 128 / 2))
    return etdrk4[2][1001][8:8:1024][ind]
end
function etdrk4_small_contour_value(x)
    ind = Int64(round(x * 128 / 2))
    return etdrk4_small_contour[2][1001][8:8:1024][ind]
end
function package_value(x)
    ind = Int64(round(x * 128 / 2))
    return package[2][21][ind]
end

using Plots
using PyPlot
pyplot()
pygui(true)
Plots.plot(x_series, [etdrk4_value, package_value, etdrk4_small_contour_value])

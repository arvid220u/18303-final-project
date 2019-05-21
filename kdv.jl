# solution of the kdv equation
# u_t = - u u_x - δ^2 u_xxx
# computation is based on v = fft(u), so linear term is diagonal

using FFTW
using Statistics
using DifferentialEquations


# solves the kdv equation using the ETDRK4 method
# returns a tuple of tt, uu with tt being timepoints and uu being the function at those timepoints
function kdv_etdrk4(h, contour_size)

    # problem specifications
    u0(x) = cos(x*pi)
    l = 2
    space_interval = (0.0,l)
    time_interval = (0.0,2)
    delta = 0.022

    # spatial discretization
    N = 1024
    x = space_interval[1] .+ space_interval[2] * (1:N)/N
    u = u0.(x)
    v = fft(u)

    # precompute ETDRK4 coefficients
    # here we utilize that things are diagonal, which means that we can store
    # the matrices as vectors instead
    k = [0:(N/2-1); 0; (-N/2+1):-1] * (2 * pi / l)# wave numbers
    L = delta ^ 2 * 1im * k .^ 3 # the linear operator for the kdv equation
    E = exp.(h*L)
    E2 = exp.(h*L/2)
    M = 32 # number of points for complex mean. need not only real part this time!
    r = contour_size * exp.(1im*2*pi*((1:M).-.5)/M); # roots of unity
    # LR: 2d array with each row representing one eigenvalue, and the different points to be evaluated for it
    LR = h*reshape(repeat(L,M), (N,M)) + reshape(repeat(r,N), (M,N))'
    Q = h*(mean( (exp.(LR/2).-1)./LR , dims = 2))
    f1 = h*(mean( (-4 .- LR .+ exp.(LR) .* (4 .- 3*LR+LR.^2))./LR.^3 , dims = 2))
    f2 = h*(mean( (2 .+LR+exp.(LR).*(-2 .+ LR))./LR.^3 , dims = 2))
    f3 = h*(mean( (-4 .-3*LR-LR.^2+exp.(LR) .* (4 .-LR))./LR.^3 , dims = 2))

    # main time-stepping loop
    uu = [u]
    tt = [time_interval[1]]
    tmax = time_interval[2]
    nmax = round(tmax/h)
    nplt = 1
    # the 0.5 comes from (u^2)_x = 2uu_x
    g = -0.5im*k
    for n in 1:nmax
        t = n*h
        # this is where we are doing the runge kutta
        Nv = g.*fft(real(ifft(v)).^2)
        a = E2.*v + Q.*Nv
        Na = g.*fft(real(ifft(a)).^2)
        b = E2.*v + Q.*Na
        Nb = g.*fft(real(ifft(b)).^2)
        c = E2.*a + Q.*(2*Nb-Nv)
        Nc = g.*fft(real(ifft(c)).^2)
        v = E.*v + Nv.*f1 + 2*(Na+Nb).*f2 + Nc.*f3

        u = reshape(real(ifft(v)), N)
        push!(uu,u)
        push!(tt,t)
    end

    return (tt, uu)
end



# solves the kdv equation using the DifferentialEquations package, after Fourier discretization
# returns a tuple (tt, uu) with tt being timepoints and uu being solutions at those timepoints
function kdv_package()

    # problem specifications
    u0(x) = cos(x*pi)
    l = 2
    space_interval = (0.0,l)
    time_interval = (0.0,2.0)
    delta = 0.022

    # spatial discretization
    N = 128
    x = space_interval[1] .+ space_interval[2] * (1:N)/N
    u = u0.(x)
    v = fft(u)

    k = [0:(N/2-1); 0; (-N/2+1):-1] * (2 * pi / l)

    function kdv_spectral(du, u, p, t)
        # geometry of u: a vector
        # first half containing u, second half containing v
        # u and v are discretized in space in lexicographical order
        # p will contain a discretized f (lexicographical), A, B, α and Δ
        k, delta, l = p
        L = delta ^ 2 * 1im * k .^ 3
        du[1:l] = L .* u - 0.5im*k .* fft(real(ifft(u)).^2)
    end

    # create the parameters
    l = N
    p = [k, delta, l]
    initial = v
    tspan = time_interval

    # create the ODEProblem
    prob = ODEProblem(kdv_spectral, initial, tspan, p)

    sol = solve(prob, saveat=0.1)

    return (sol.t, real.(ifft.(sol.u)))
end

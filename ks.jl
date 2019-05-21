
using FFTW
using Statistics

# returns a tuple (tt,uu) containing the timesteps and the functionvalues at those timesteps
function ks_etdrk4(h)
    # solution of Kuramoto-Sivashinsky
    # u_t = -u*u_x - u_xx - u_xxxx, periodic BCs on [0,32*pi]
    # computation is based on v = fft(u), so linear term is diagonal
    # the code below is inspired by Kassam and Trefethen's Matlab code

    # problem specifications
    u0(x) = cos(x/16) * (1 + sin(x/16))
    l = 32*pi
    space_interval = (0.0,l)
    time_interval = (0.0,150.0)

    # spatial discretization
    N = 128
    x = space_interval[1] .+ space_interval[2] * (1:N)/N
    u = u0.(x)
    v = fft(u)

    # precompute ETDRK4 coefficients
    # here we utilize that things are diagonal, which means that we can store
    # the matrices as vectors instead
    k = [0:(N/2-1); 0; (-N/2+1):-1] * (2 * pi) / l # wave numbers
    L = k .^ 2 - k .^ 4
    E = exp.(h*L)
    E2 = exp.(h*L/2)
    M = 16 # number of points for complex mean
    r = exp.(1im*pi*((1:M).-.5)/M); # roots of unity
    # LR: 2d array with each row representing one eigenvalue, and the different points to be evaluated for it
    LR = h*reshape(repeat(L,M), (N,M)) + reshape(repeat(r,N), (M,N))'
    Q = h*real(mean( (exp.(LR/2).-1)./LR , dims = 2))
    f1 = h*real(mean( (-4 .- LR .+ exp.(LR) .* (4 .- 3*LR+LR.^2))./LR.^3 , dims = 2))
    f2 = h*real(mean( (2 .+LR+exp.(LR).*(-2 .+ LR))./LR.^3 , dims = 2))
    f3 = h*real(mean( (-4 .-3*LR-LR.^2+exp.(LR) .* (4 .-LR))./LR.^3 , dims = 2))

    # main time-stepping loop
    uu = [u]
    tt = [time_interval[1]]
    tmax = time_interval[2]
    nmax = round(tmax/h)
    nplt = floor((tmax/100)/h);
    g = -0.5im*k
    for n in 1:nmax
        t = n*h
        Nv = g.*fft(real(ifft(v)).^2)
        a = E2.*v + Q.*Nv
        Na = g.*fft(real(ifft(a)).^2)
        b = E2.*v + Q.*Na
        Nb = g.*fft(real(ifft(b)).^2)
        c = E2.*a + Q.*(2*Nb-Nv)
        Nc = g.*fft(real(ifft(c)).^2)
        v = E.*v + Nv.*f1 + 2*(Na+Nb).*f2 + Nc.*f3
        if mod(n,nplt) == 0
            u = reshape(real(ifft(v)), N)
            push!(uu,u)
            push!(tt,t)
        end
    end

    return (tt, uu)
end

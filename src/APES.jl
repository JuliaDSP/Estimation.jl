@doc """
Amplitude and Phase Estimation of a Sinusoid (APES).  
Implemented as described in [Li, Li & Stoica 1998](http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=700967).

This function estimates the complex amplitude of signal `y` at frequencies `f`.  
M denotes the filter length, which is set to N/2 by default [[Li & Stoica 1996](http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=506612)].

## References
H. Li, J. Li, and P. Stoica. Performance analysis of forward-backward matched-filterbank spectral esti- mators. Signal Processing, IEEE Transactions on, 46(7):1954–1966, 1998.

J. Li and P. Stoica. An adaptive filtering approach to spectral estimation and sar imaging. Signal Processing, IEEE Transactions on, 44(6):1469–1484, 1996.

Y. Hua, A. Gershman, and Q. Cheng. High-resolution and robust signal processing. Signal Processing and Communications. CRC Press, October 2003.

""" ->
function apes{T <: FloatingPoint}(y::Array{T}, f::Union(AbstractArray{T}, T); M::Int=div(length(y), 2), Q::Union(Array{T, 2}, Nothing)=nothing)

    N = length(y)
    L = N - M + 1  # (eqn 3)

    y_f = zeros(M, L)  # Forward data vectors
    y_b = zeros(M, L)  # Backward data vectors
    R_f = zeros(M, M)  # Forward data covariance matrix
    R_b = zeros(M, M)  # Backward data covariance matrix
    for l = 0:L-1
        y_f[:, l+1] = y[l+1 : 1 : l+M-1+1]               # eqn 3
        y_b[:, l+1] = conj(y[N-l-1+1 : -1 : N-l-M+1])    # eqn 6
        R_f += y_f[:, l+1] * y_f[:, l+1]'                # eqn 12
        R_b += y_b[:, l+1] * y_b[:, l+1]'                # eqn 13
    end
    R_f = R_f ./ L
    R_b = R_b ./ L

    # Forward-backward estimate of covariance matrix (eqn 14)
    R  = (1/2) * (R_f + R_b)

    # Calculate complex amplitude for each frequency
    α = zeros(Complex{T}, size(f))
    for i in 1:length(f)

        ω = 2 * pi * f[i]
        α[i] = apes(ω, M, L, y_f, y_b, R, Q)
    end

    return α
end

function apes{T <: FloatingPoint}(ω::T, M::Int, L::Int, y_f::Array{T, 2}, y_b::Array{T, 2}, R::Array{T, 2}, Q::Union(Nothing, Array{T, 2}))

    a_m = exp(-im * ω .* [0:1:M-1 ; ])
    a_l = exp(-im * ω .* [0:1:L-1 ; ])

    # Normalised Fourier transforms of forward and backward data vectors
    g_f = (1/L) * y_f * conj(a_l)  # eqn 4.3.11 [1]
    g_b = (1/L) * y_b * conj(a_l)
    G   = (1/sqrt(2)) * [g_f g_b]  # eqn 32

    # Estimate of noise covariance matrix
    if Q == nothing
        Q  = R - G * ctranspose(G)  # eqn 31
    end
    Qi = inv(Q)

    # Complex amplitude of sinusoid signal with frequency ω (eqn 34)
    (2 * (a_m' * Qi * g_f) / (a_m' * Qi * a_m))[1]
end

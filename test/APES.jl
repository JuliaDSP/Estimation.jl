#
# Test signal
#

## Li & Stoica fig1
spectral_lines_f = [0.04, 0.05, 0.08, 0.09,    0.16:0.005:0.35]
spectral_lines_p = vec(repmat([pi/4], length(spectral_lines_f)))
spectral_lines_a = [10, 10, 10, 2,   vec(repmat([0.5], length(spectral_lines_f)-4, 1))]

## Simulation properties
N = 128
noise_mean = 0
noise_vari = 1

## Create example signal
function create_example(N, f, a, p)
    z = vec(zeros(N, 1))
    n = 1.0:N
    for i in 1:length(f)
        z = z + a[i] * sin(2 * pi * f[i] * n + p[i])
    end
    return z
end

# For passing in known cov matrix
function covariance_matrix{T <: FloatingPoint}(s::Array{T}; M::Int=int(length(s)/2))
    N = length(s)
    L = N - M + 1
    CovMatrix = zeros(M, M)
    for l = 0:N-M
        CovMatrix += s[l+1 : 1 : l+M-1+1] * s[l+1 : 1 : l+M-1+1]'
    end
    CovMatrix ./ L
end

# For reproducible results
srand(1)


#
# Test example from Li & Stoica
#

signal = create_example(N, spectral_lines_f, spectral_lines_a, spectral_lines_p)
noise  = sqrt(noise_vari) * randn(size(signal))
sn     = signal + noise
MSE    = mean((abs(apes(sn, spectral_lines_f)) .- spectral_lines_a).^2)

@test_approx_eq abs(apes(sn, spectral_lines_f)[1:4])  [10.150160442095107, 9.871766839666593, 10.02406973575788, 2.131756870868149]
@test_approx_eq MSE  0.10794065421715185


#
# Test at lower SNR
#

signal  = create_example(N, spectral_lines_f, spectral_lines_a, spectral_lines_p)
noise   = sqrt(noise_vari) * randn(size(signal))
sn      = signal + noise * 50
MSE     = mean((abs(apes(sn, spectral_lines_f)) .- spectral_lines_a).^2)
MSE_w_Q = mean((abs(apes(sn, spectral_lines_f, Q = covariance_matrix(noise))) .- spectral_lines_a).^2)

@test_approx_eq MSE     92.84644647200092
@test_approx_eq MSE_w_Q 9.936179958764319

using Estimation
using Base.Test

function noiseless_bias()
	n = 65536
	N = 100

	result = []
	for trial=1:N
		freq = 0.1 + 0.3*rand()
		sig = [cos(2.0*pi*freq*i) for i=1:n]

		est_freq = Estimation.quinn_fernandes(sig, 2.0*pi*0.1)/(2.0*pi)
		result = [(freq-est_freq),result]
	end

	result
end

function noisy_bias()
	n = 65536
	N = 100

	result = []
	for trial=1:N
		freq = 0.1 + 0.3*rand()
		sig = [cos(2.0*pi*freq*i) for i=1:n] + 0.1*randn(n)

		est_freq = Estimation.quinn_fernandes(sig, 2.0*pi*0.1)/(2.0*pi)
		result = [(freq-est_freq),result]
	end

	result
end
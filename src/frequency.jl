function quinn_fernandes(signal, seed)
	eps = 1.0e-6

	n = length(signal)
	a = 2.0*cos(seed)
	b = 0.1
	signal = vcat([0,0], signal)
	xi = zeros(size(signal))

	it = 0
	while (abs(a-b) > eps)
		it = it+1
		if (it > 1000)
			b = NaN
			break
		end	
		a = 2*b - a
		for i=3:(n-1)
			xi[i] = signal[i] + a*xi[i-1] - xi[i-2]
		end
		b = sum([(xi[i]+xi[i-2])*xi[i-1] for i=3:(n-1)])/sum([xi[i-1]^2 for i=3:(n-1)])
	end

	acos(0.5*b)
end
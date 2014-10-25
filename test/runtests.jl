using Estimation
using Base.Test

sig = [cos(2*pi*0.15*i) for i=1:2048];
@test_approx_eq_eps Estimation.quinn_fernandes(sig,2*pi*0.1)/(2*pi) 0.15 1e-5
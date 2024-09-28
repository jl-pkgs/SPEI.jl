using Test

p = [1.30, 17.0]
D = Gamma(p[1], p[2])

@test cdfgam(1.0, p) == cdf(D, 1.0)

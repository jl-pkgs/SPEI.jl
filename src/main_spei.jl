function spei(x::AbstractVector)
  x2 = x[.!isnan(x)] |> sort
  beta = pwm(x2, 0.0, 0.0, 0)
  params = _fit_logLogistic(beta)

  n = length(x)
  z = zeros(Float64, n)
  for i in 1:n
    prob = cdf_logLogistic(x[i], params)
    z[i] = -invcdf_standardGaussian(prob)
  end
  (; z, coef=params)
end


# 这里需要提供一个reference period
"""
    spi(x::AbstractVector; fit="lmom")

# Arguements
- `fit`: 
  + `lmom`: use lmoments to fit gamma distribution, R version
  + `mle`: maximum likelihood estimation, julia solver
"""
function spi(x::AbstractVector; fit="lmom")
  x2 = x[.!isnan(x)]
  x_nozero = x2[x2.>0]
  q = 1 - length(x_nozero) / length(x2) # probability of zero
  
  if fit == "mle"
    D = fit_mle(Gamma, x_nozero)
  else
    p = fit_gamma(x_nozero)
    D = Gamma(p...)
  end

  z = fill(NaN, length(x))
  for i in eachindex(z)
    if !isnan(x[i])
      p = q .+ (1 - q) * cdf(D, x[i])
      z[i] = qnorm(p)
      x[i] < 0 && (z[i] = -Inf)
    end
  end
  (; z, coef=params(D))
end

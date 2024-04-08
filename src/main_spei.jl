function param_spei(x_ref::AbstractVector)
  x2 = x_ref[.!isnan.(x_ref)] |> sort # have to sort
  beta = pwm(x2, 0.0, 0.0, 0)
  params = _fit_logLogistic(beta)
  params
end


# add a reference function to spi
function param_spi(x_ref::AbstractVector; fit="lmom")
  x2 = x_ref[.!isnan.(x_ref)]
  x_nozero = x2[x2.>0]
  q = 1 - length(x_nozero) / length(x2) # probability of zero

  if fit == "mle"
    D = fit_mle(Gamma, x_nozero)
  else
    p = fit_gamma(x_nozero)
    D = Gamma(p...)
  end
  params(D), q
end


"""
    spei(x::AbstractVector)
"""
function spei(x::AbstractVector, x_ref::AbstractVector=x)
  params = param_spei(x_ref)
  z = zeros(Float64, length(x))

  for i in eachindex(z)
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
function spi(x::AbstractVector, x_ref::AbstractVector=x; fit="lmom")
  param, q = param_spi(x_ref; fit)
  D = Gamma(param...)
  
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

# TODO: 需要添加scale参数

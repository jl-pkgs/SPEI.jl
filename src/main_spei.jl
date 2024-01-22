function spei(x::AbstractVector)
  x2 = x[.!isnan(x)] |> sort
  beta = pwm(x2, 0.0, 0.0, 0)
  params = fit_logLogistic(beta)

  n = length(x)
  z = zeros(Float64, n)
  for i in 1:n
    prob = cdf_logLogistic(x[i], params)
    z[i] = -invcdf_standardGaussian(prob)
  end
  (; z, coef=params)
end


# 这里需要提供一个reference period
function spi(x::AbstractVector)
  x2 = x[.!isnan(x)]
  x_nozero = x2[x2.>0]
  q = 1 - length(x_nozero) / length(x2) # probability of zero
  D = fit_mle(Gamma, x_nozero)

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


function spi_c(x::AbstractVector)
  x2 = x[.!isnan(x)]
  x_nozero = x2[x2.>0]
  q = 1 - length(x_nozero) / length(x2) # probability of zero

  beta = pwm(x_nozero, 0.0, 0.0, 0)
  
  params = fit_gamma(beta)
  z = fill(NaN, length(x))
  for i = eachindex(z)
    if x[i] > 0
      z[i] = q_gamma(x[i], params)
    end
  end

  (; z, coef=params)
end

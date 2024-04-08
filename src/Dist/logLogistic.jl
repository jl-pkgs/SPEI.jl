function _fit_logLogistic(beta::AbstractVector)
  # params = zeros(3)
  g1 = 0.0
  g2 = 0.0

  # estimate gamma parameter
  γ = (2 * beta[2] - beta[1]) / (6 * beta[2] - beta[1] - 6 * beta[3])
  g1 = exp(loggamma(1 + 1 / γ))
  g2 = exp(loggamma(1 - 1 / γ))

  # estimate alpha parameter
  α = (beta[1] - 2 * beta[2]) * γ / (g1 * g2) # params[2]
  # estimate beta parameter
  β = beta[1] - α * g1 * g2                   # params[1]
  β, α, γ
end


function fit_logLogistic(x::AbstractVector)
  x2 = x[.!isnan.(x)] |> sort
  beta = pwm(x2, 0.0, 0.0, 0)
  params = _fit_logLogistic(beta)
  params
end


function cdf_logLogistic(x::Real, params)
  β, α, γ = params[1:3]
  1 / (1 + (pow(α / (x - β), γ)))
end


function invcdf_standardGaussian(prob::Float64)
  C = [2.515517, 0.802853, 0.010328]
  d = [0, 1.432788, 0.189269, 0.001308]

  W = prob <= 0.5 ? sqrt(-2 * log(prob)) : sqrt(-2 * log(1 - prob))
  WW = W * W
  WWW = WW * W
  resul = W - (C[1] + C[2] * W + C[3] * WW) / (1 + d[2] * W + d[3] * WW + d[4] * WWW)

  prob > 0.5 && (resul = -resul)
  return resul
end

export fit_logLogistic

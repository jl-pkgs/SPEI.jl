function fit_logLogistic(beta::AbstractVector)
  params = zeros(3)
  g1 = 0.0
  g2 = 0.0

  # estimate gamma parameter
  params[3] = (2 * beta[2] - beta[1]) / (6 * beta[2] - beta[1] - 6 * beta[3])
  g1 = exp(loggamma(1 + 1 / params[3]))
  g2 = exp(loggamma(1 - 1 / params[3]))

  # estimate alpha parameter
  params[2] = (beta[1] - 2 * beta[2]) * params[3] / (g1 * g2)
  # estimate beta parameter
  params[1] = beta[1] - params[2] * g1 * g2
  params
end


function cdf_logLogistic(value::Float64, params::AbstractVector)
  return 1 / (1 + ((params[2] / (value - params[1]))^params[3]))
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

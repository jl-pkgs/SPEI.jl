function _fit_logLogistic(β1::T, β2::T, β3::T) where {T<:Real}
  g1 = 0.0
  g2 = 0.0

  # estimate gamma parameter
  γ = (2 * β2 - β1) / (6 * β2 - β1 - 6 * β3)
  g1 = exp(loggamma(1 + 1 / γ))
  g2 = exp(loggamma(1 - 1 / γ))

  # estimate alpha parameter
  α = (β1 - 2 * β2) * γ / (g1 * g2) # params[2]
  # estimate beta parameter
  β = β1 - α * g1 * g2                   # params[1]
  β, α, γ
end


function fit_logLogistic(x::AbstractVector)
  x2 = x[.!isnan.(x)] |> sort
  β1, β2, β3 = pwm(x2, 0.0, 0.0, 0)
  params = _fit_logLogistic(β1, β2, β3)
  params
end


function cdf_logLogistic(x::Real, params)
  β, α, γ = params[1:3]
  1 / (1 + (pow(α / (x - β), γ)))
end


const _C = [2.515517, 0.802853, 0.010328]
const _d = [0.0, 1.432788, 0.189269, 0.001308]

function invcdf_standardGaussian(prob::T) where {T<:Real}
  W = prob <= 0.5 ? sqrt(-2 * log(prob)) : sqrt(-2 * log(1 - prob))
  WW = W * W
  WWW = WW * W
  resul = W - (_C[1] + _C[2] * W + _C[3] * WW) / (1 + _d[2] * W + _d[3] * WW + _d[4] * WWW)

  prob > 0.5 && (resul = -resul)
  return resul
end

export fit_logLogistic

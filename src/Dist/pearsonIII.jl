function pearsonIIIFit(L::Array{Float64,1})
  params = zeros(3)
  sqrtPi = sqrt(pi)
  lmom_ratio = L[4] / L[3]  # third l-Moment ratio

  if lmom_ratio > 0.0 && lmom_ratio < 0.333333
    t = 3 * pi * lmom_ratio * lmom_ratio
    t2 = t * t
    t3 = t2 * t
    params[3] = (1 + 0.2906 * t) / (t + 0.1882 * t2 + 0.0442 * t3)
  elseif lmom_ratio >= 0.333333 && lmom_ratio < 1
    t = 1 - lmom_ratio
    t2 = t * t
    t3 = t2 * t
    params[3] = (0.36067 * t - 0.5967 * t2 + 0.25361 * t3) /
                          (1 - 2.78861 * t + 2.56096 * t2 - 0.77045 * t3)
  end
  params[2] = sqrtPi * L[2] *
                        exp(loggamma(params[3]) - loggamma(params[3] + 0.5))
  params[1] = L[1] - params[2] * params[3]
  params
end

function pearsonIIIStandardize(value::Float64, params::Array{Float64,1})
  origin = params[1]  # location
  alfa = params[2]    # shape
  beta = params[3]    # scale
  z = 0.0

  if value >= origin
    y = (value - origin) / alfa
    z = (pow(y / beta, 0.333333333333333) + 1 / (9 * beta) - 1) * pow(9 * beta, 0.5)
  else
    z = -3.09023230616779
  end
  
  return z
end

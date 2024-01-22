function gammaFit(L::AbstractVector)
  params = zeros(2)
  lmom_ratio = L[3] / L[2]  # second l-Moment ratio
  
  if lmom_ratio > 0.0 && lmom_ratio < 0.5
    t = pi * lmom_ratio * lmom_ratio
    t2 = t * t
    t3 = t2 * t
    params[2] = (1 - 0.308 * t) / (t - 0.05812 * t2 + 0.01765 * t3)
    params[1] = L[1] / params[2]
  elseif lmom_ratio >= 0.5 && lmom_ratio < 1
    t = 1 - lmom_ratio
    t2 = t * t
    params[2] = (0.7213 * t - 0.5947 * t2) / (1 - 2.1817 * t + 1.2113 * t2)
    params[1] = L[2] / params[2]
  end
  params
end


function q_gamma(value::Float64, params::AbstractVector)
  alfa = params[1]
  beta = params[2]
  
  y = value / alfa / beta
  z = (y^(1 / 3) + 1 / (9 * beta) - 1) * sqrt(9 * beta)
  return z
end

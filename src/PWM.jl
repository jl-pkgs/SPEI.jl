# https://github.com/cran/TLMoments/blob/master/src/PWM.cpp
function PWM(x::AbstractVector, r::Int)
  n = length(x)
  xs = sort(x)
  w = 0.0
  vorfaktor = 1.0

  for k = n:-1:(n-r)
    vorfaktor /= k
  end

  sum = 0.0
  for j = 1:n
    w = 1.0
    for i = 1:r
      w *= (j - i)
    end
    sum += w * xs[j]
  end
  return sum * vorfaktor
end

PWM(x::AbstractVector, r::AbstractVector{Int}) = [PWM(x, ri) for ri in r]


function pwm2lmom(pwm::AbstractVector)
  nmom = length(pwm)
  L = zeros(nmom)
  R = zeros(nmom)
  R[1] = NaN

  for i in 1:nmom
    r = i - 1
    sum = 0
    for k in 0:r
      weight = (-1)^(r - k) * binomial(r, k) * binomial(r + k, k)
      sum += weight * pwm[k+1]
    end
    L[i] = sum
  end

  nmom >= 2 && (R[2] = L[2] / L[1])
  for r in 3:nmom
    R[r] = L[r] / L[2]
  end
  (; lambdas=L, ratios=R, source="pwm2lmom")
end

export PWM, pwm2lmom

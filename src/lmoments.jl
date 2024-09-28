"""
Calculates the first three probability weighted moments of a sample, using
either the unbiased estimator (when A=B=0) or a plotting position formula (when
A<=B<=0).

This are alpha PWMs, following Rao & Hamed 2000, eqs. 3.1.4 and 3.1.6
"""
function pwm(series::AbstractVector{T}, A::Float64, B::Float64, isBeta::Int) where {T<:Real}
  F = 0.0
  n = length(series)
  acum1 = acum2 = acum3 = T(0)

  if A == 0 && B == 0  # use unbiased estimator
    for i = 1:n
      acum1 += series[i]
      if isBeta == 0  # compute alpha PWMs
        acum2 += series[i] * (n - i)
        acum3 += series[i] * (n - i) * (n - i - 1)
      end
      if isBeta == 1  # compute beta PWMs
        acum2 += series[i] * (i - 1)
        acum3 += series[i] * (i - 1) * (i - 2)
      end
    end
  else  # use plotting-position (biased) estimator
    for i = 1:n
      acum1 += series[i]
      F = (i + A) / (n + B)
      if isBeta == 0  # compute alpha PWMs
        acum2 += series[i] * (1 - F)
        acum3 += series[i] * (1 - F) * (1 - F)
      end
      if isBeta == 1  # compute beta PWMs
        acum2 += series[i] * F
        acum3 += series[i] * F * F
      end
    end
  end

  β1 = acum1 / n
  β2 = acum2 / n / (n - 1)
  β3 = acum3 / n / ((n - 1) * (n - 2))
  β1, β2, β3
end

function lmoments(series::AbstractVector{T}, A::Float64, B::Float64) where {T<:Real}
  # alpha = zeros(Float64, 3)
  # Calculate the first three PWMs
  lmom = zeros(T, 3)
  alpha = pwm(series, A, B, 0)

  lmom[1] = alpha[1]
  lmom[2] = alpha[1] - 2 * alpha[2]
  lmom[3] = alpha[1] - 6 * alpha[2] + 6 * alpha[3]
  # lMoment[4] = alpha[1] - 12*alpha[2] + 30*alpha[3] - 20*alpha[4]
  lmom
end

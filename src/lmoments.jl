"""
Calculates the first three probability weighted moments of a sample, using
either the unbiased estimator (when A=B=0) or a plotting position formula (when
A<=B<=0).

This are alpha PWMs, following Rao & Hamed 2000, eqs. 3.1.4 and 3.1.6
"""
function pwm!(acum::AbstractVector{T}, pwms::AbstractVector{T},
  series::AbstractVector{T}, A::Float64, B::Float64, isBeta::Int) where {T<:Real}

  F = 0.0
  n = length(series)

  if A == 0 && B == 0  # use unbiased estimator
    for i = 1:n
      acum[1] += series[i]
      if isBeta == 0  # compute alpha PWMs
        acum[2] += series[i] * (n - i)
        acum[3] += series[i] * (n - i) * (n - i - 1)
      end
      if isBeta == 1  # compute beta PWMs
        acum[2] += series[i] * (i - 1)
        acum[3] += series[i] * (i - 1) * (i - 2)
      end
    end
  else  # use plotting-position (biased) estimator
    for i = 1:n
      acum[1] += series[i]
      F = (i + A) / (n + B)
      if isBeta == 0  # compute alpha PWMs
        acum[2] += series[i] * (1 - F)
        acum[3] += series[i] * (1 - F) * (1 - F)
      end
      if isBeta == 1  # compute beta PWMs
        acum[2] += series[i] * F
        acum[3] += series[i] * F * F
      end
    end
  end

  pwms[1] = acum[1] / n
  pwms[2] = acum[2] / n / (n - 1)
  pwms[3] = acum[3] / n / ((n - 1) * (n - 2))
  pwms
end

function pwm(series::AbstractVector{T}, A::Float64, B::Float64, isBeta::Int) where {T<:Real}
  acum = zeros(T, 3)
  pwms = zeros(T, 3)
  pwm!(acum, pwms, series, A, B, isBeta)
end

function lmoments(series::AbstractVector{T}, A::Float64, B::Float64) where { T<:Real }
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

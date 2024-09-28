function param_spei(x_ref::AbstractVector)
  lgl = .!isnan.(x_ref)
  x2 = x_ref[lgl]
  sort!(x2) # have to sort

  β1, β2, β3 = pwm(x2, 0.0, 0.0, 0) # improve at here
  params = _fit_logLogistic(β1, β2, β3)
  params
end

# add a reference function to spi
function param_spi(x_ref::AbstractVector; fit="lmom")
  lgl = .!isnan.(x_ref)
  x2 = @view x_ref[lgl]
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
    spei(x::AbstractVector, x_ref::AbstractVector=x)
    spei!(z::AbstractVector{T}, x::AbstractVector, x_ref::AbstractVector=x; nmin=3) where {T<:Real}
"""
function spei!(z::AbstractVector{T}, x::AbstractVector, x_ref::AbstractVector=x; nmin=3) where {T<:Real}
  # 如果有效观测数量不够，需要跳过
  n_valid = 0
  for xi in x_ref
    xi == xi && (n_valid += 1)
  end
  n_valid < nmin && (return T(NaN), T(NaN), T(NaN)) # β, α, γ

  params = param_spei(x_ref)
  @inbounds for i in eachindex(z)
    prob = cdf_logLogistic(x[i], params)
    z[i] = -invcdf_standardGaussian(prob)
  end
  params
end


"""
    spi(x::AbstractVector; fit="lmom")

# Arguements
- `fit`: 
  + `lmom`: use lmoments to fit gamma distribution, R version
  + `mle`: maximum likelihood estimation, julia solver
"""
function spi!(z::AbstractVector, x::AbstractVector, x_ref::AbstractVector=x; nmin=3, fit="lmom")
  n_valid = 0
  for xi in x_ref
    xi > 0 && (n_valid += 1)
  end
  n_valid < nmin && (return T(NaN), T(NaN)) # β, α, γ

  param, q = param_spi(x_ref; fit)
  D = Gamma(param...)

  @inbounds for i in eachindex(z)
    if !isnan(x[i])
      p = q .+ (1 - q) * cdf(D, x[i])
      z[i] = qnorm(p)
      x[i] < 0 && (z[i] = -Inf)
    end
  end
  params(D)
end


function spei(x::AbstractVector{T}, x_ref::AbstractVector{T}=x; nmin::Int=3) where {T<:Real}
  z = fill(T(NaN), length(x))
  params = spei!(z, x, x_ref; nmin)
  (; z, coef=params)
end

function spi(x::AbstractVector{T}, x_ref::AbstractVector{T}=x; nmin::Int=3, fit="lmom") where {T<:Real}
  z = fill(T(NaN), length(x))
  params = spi!(z, x, x_ref; nmin, fit)
  (; z, coef=params)
end


export spei!, spi!, spei, spi

using SpecialFunctions: gamma_inc, lgamma
export cdfgam, pelgam


function cdfgam(x::Float64, para::AbstractVector{Float64})
  α, β = para
  if α <= 0 || β <= 0
    println(" *** ERROR *** ROUTINE CDFGAM : PARAMETERS INVALID")
    return NaN
  end
  x <= 0 && return NaN
  gamma_inc(α, x / β)[1]
end


function pelgam(xmom::AbstractVector{Float64})
  a1, a2, a3 = -0.3080, -0.05812, 0.01765
  b1, b2, b3, b4 = 0.7213, -0.5947, -2.1817, 1.2113

  if xmom[1] <= xmom[2] || xmom[2] <= 0.0
    println(" *** ERROR *** ROUTINE PELGAM : L-MOMENTS INVALID")
    return 0.0, 0.0
  end

  cv = xmom[2] / xmom[1]
  if cv < 0.5
    t = pi * cv * cv
    alpha = (1.0 + a1 * t) / (t * (1.0 + t * (a2 + t * a3)))
  else
    t = 1.0 - cv
    alpha = t * (b1 + t * b2) / (1.0 + t * (b3 + t * b4))
  end
  alpha, xmom[1] / alpha # α, β
end


function lmrgam(para::Vector{Float64}, nmom::Int=2)
  CONST = 0.564189583547756287
  a0 = 0.32573501
  a1, a2, a3 = 0.16869150, 0.78327243, -0.29120539
  b1, b2 = 0.46697102, 0.24255406
  c0 = 0.12260172
  c1, c2, c3 = 0.53730130, 0.43384378, 0.11101277
  d1, d2 = 0.18324466, 0.20166036
  e1, e2, e3 = 0.23807576, 0.15931792, 0.11618371
  f1, f2, f3 = 0.51533299, 0.71425260, 0.19745056
  g1, g2, g3 = 0.21235833, 0.41670213, 0.31925299
  h1, h2, h3 = 0.90551443, 0.26649995, 0.26193668

  alpha, beta = para[1:2]
  if alpha <= 0.0 || beta <= 0.0
    println(" *** ERROR *** ROUTINE LMRGAM : PARAMETERS INVALID")
    return
  end

  if nmom > 4
    println(" *** ERROR *** ROUTINE LMRGAM : PARAMETER NMOM TOO LARGE")
    return
  end

  xmom = zeros(nmom)
  xmom[1] = alpha * beta
  nmom == 1 && return xmom

  xmom[2] = beta * CONST * exp(lgamma(alpha + 0.5) - lgamma(alpha))
  nmom == 2 && return xmom

  if alpha >= 1.0
    z = 1.0 / alpha
    xmom[3] = sqrt(z) * (((a3 * z + a2) * z + a1) * z + a0) / ((b2 * z + b1) * z + 1.0)
    nmom == 3 && return xmom
    
    xmom[4] = (((c3 * z + c2) * z + c1) * z + c0) / ((d2 * z + d1) * z + 1.0)
  else
    z = alpha
    xmom[3] = (((e3 * z + e2) * z + e1) * z + 1.0) / (((f3 * z + f2) * z + f1) * z + 1.0)
    nmom == 3 && return xmom

    xmom[4] = (((g3 * z + g2) * z + g1) * z + 1.0) / (((h3 * z + h2) * z + h1) * z + 1.0)
  end

  return xmom
end

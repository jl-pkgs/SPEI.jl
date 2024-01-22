module SPEI

using Distributions

export spei, spi

include("math.jl")
include("lmoments.jl")

include("Dist/Dist.jl")


function spei(x::AbstractVector)
  x2 = x[.!isnan(x)] |> sort
  beta = pwm(x2, 0.0, 0.0, 0)
  params = logLogisticFit(beta)

  n = length(x)
  z = zeros(Float64, n)
  for i in 1:n
    prob = logLogisticCDF(x[i], params)
    z[i] = -standardGaussianInvCDF(prob)
  end
  (; z, coef=params)
end


function spi(x::AbstractVector)
  x2 = x[.!isnan(x)]
  x_nozero = x2[x2 .> 0] 
  q = 1 - length(x_nozero) / length(x2) # probability of zero
  D = fit_mle(Gamma, x_nozero)
  
  z = fill(NaN, length(x))  
  for i in eachindex(z)
    if !isnan(x[i])
      p = q .+ (1 - q) * cdf(D, x[i])
      z[i] = quantile(Normal(), p)
      x[i] < 0 && (z[i] = -Inf)
    end
  end
  
  (; z, D=D)
end


end # module SPEI

module SPEI

export spei, spi, spi_c
export param_spei, param_spi
export qnorm

using Distributions
using SpecialFunctions: gamma_inc, lgamma, loggamma


qnorm(p::Real) = quantile(Normal(), p)

export pow
# pow = ^
# import Base:^
pow(x, y) = x^y
function pow(x::T, y::AbstractFloat) where {T<:AbstractFloat}
  x < 0 ? T(NaN) : x^y
end


include("PWM.jl")
include("lmoments.jl")

include("DIST/DIST.jl")
include("main_spei.jl")

end # module SPEI

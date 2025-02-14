module SPEI

export spei, spi, spi_c
export param_spei, param_spi
export qnorm

using Distributions
using SpecialFunctions: gamma_inc, lgamma, loggamma
using Ipaper, Dates
import NaNStatistics:nanmean, nanstd

function get_dn(date; delta=8)
  days = Dates.dayofyear(date)
  cld(days, delta) # int
end

qnorm(p::Real) = quantile(Normal(), p)

export pow
# pow = ^
# import Base:^
# "Faster method for exponentiation"
# @fastmath pow(x::Real, y::Real) = exp(y * log(x))
@fastmath pow(x, y) = exp(y * log(x))

@fastmath function pow(x::T, y::AbstractFloat) where {T<:AbstractFloat}
  x < 0 ? T(NaN) : exp(y * log(x))
end

include("PWM.jl")
include("lmoments.jl")

include("DIST/gamma.jl")
include("DIST/logLogistic.jl")
# include("DIST/DIST.jl")
include("main_spei.jl")

include("drought_ZSI.jl")
include("drought_SPEI.jl")

end # module SPEI

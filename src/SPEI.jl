module SPEI

export spei, spi, spi_c
export qnorm

using Distributions
using SpecialFunctions: gamma_inc, lgamma, loggamma


qnorm(p::Real) = quantile(Normal(), p)
pow = ^


include("PWM.jl")
include("lmoments.jl")

include("Dist/Dist.jl")
include("main_spei.jl")

end # module SPEI

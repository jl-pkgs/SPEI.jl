module SPEI

export spei, spi, spi_c

using Distributions
using SpecialFunctions: gamma_inc, lgamma, loggamma


include("math.jl")
include("PWM.jl")
include("lmoments.jl")

include("Dist/Dist.jl")
include("main_spei.jl")

end # module SPEI

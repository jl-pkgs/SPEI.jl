module SPEI

using Distributions

export spei, spi, spi_c

include("math.jl")
include("lmoments.jl")

include("Dist/Dist.jl")
include("main_spei.jl")

include("Lmoments/gamma.jl")

end # module SPEI

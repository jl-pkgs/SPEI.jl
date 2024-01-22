module SPEI

export spei

include("math.jl")
include("lmoments.jl")
include("logLogistic.jl")
include("pearsonIII.jl")
include("gamma.jl")

function spei(x; scale=6)
  sort!(x)
  beta = pwm(x, 0.0, 0.0, 0)
  params = logLogisticFit(beta)
  prob = logLogisticCDF(x[1], params)
  z = -standardGaussianInvCDF(prob)
  z
end


end # module SPEI

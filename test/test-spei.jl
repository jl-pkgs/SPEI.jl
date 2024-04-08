using Test
using RCall
using SPEI

R"""
library(SPEI2)
"""

function nanmax(x)
  x2 = x[.!isnan.(x)]
  maximum(x2)
end

# wb = R"wb" |> rcopy |> x -> Float32.(x)
wb = deserialize("data/wb")
x, x_ref = deserialize("data/x_ref")

@testset "spei" begin
  @time r_R = R"cal_spei($x, $x_ref)$z" |> rcopy
  @time r_jl = spei(x, x_ref).z

  # z_r = r_R[:z]
  # z_jl = r_jl.z
  e = r_R - r_jl
  @test nanmax(e) <= 1e-3
end

@testset "spi" begin
  @time r_R = R"cal_spi(wb)" |> rcopy
  @time r_jl = spi(wb)

  z_r = r_R[:z]
  z_jl = r_jl.z
  e = z_r - z_jl
  @test abs(nanmax(e)) <= 1e-6
end

# > 计算结果完全一致
# begin
#   using Plots
#   plot(z_r, label="R")
#   plot!(z_jl, label="Julia")
#   plot!(e, label="R - Julia")
# end

using Test
using RCall


R"""
library(SPEI2)
"""
wb = R"wb" |> rcopy


@testset "spei" begin
  @time r_R = R"cal_spei(wb)" |> rcopy
  @time r_jl = spei(wb)

  z_r = r_R[:z]
  z_jl = r_jl.z
  e = z_r - z_jl
  @test maximum(abs.(e)) <= 1e-3
end


@testset "spi" begin
  @time r_R = R"cal_spi(wb)" |> rcopy
  @time r_jl = spi(wb)

  z_r = r_R[:z]
  z_jl = r_jl.z
  e = z_r - z_jl

end

## 误差还挺大的
begin
  using Plots

  plot(z_r, label="R")
  plot!(z_jl, label="Julia")
  plot!(e, label="R - Julia")
end

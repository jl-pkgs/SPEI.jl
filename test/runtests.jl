using Test
using SPEI
using Serialization
nanmaximum(x) = maximum(x[.!isnan.(x)])

dir_proj = "$(@__DIR__)/.."
wb = deserialize("$dir_proj/data/wb")
x = wb[wb.>0]

include("test-drought.jl")

@testset "gamma and loglogistic" begin
  @test fit_gamma(x) == (1.2729799501897179, 17.435758906451866)
  @test fit_logLogistic(x) == (-22.35172530856679, 40.88689577169438, 4.418130441097305)
end

@testset "spei" begin
  @test spei(x).coef == (-22.35172530856679, 40.88689577169438, 4.418130441097305)
  @test spi(x).coef == (1.2729799501897179, 17.435758906451866)
end

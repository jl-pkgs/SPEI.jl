using Test
using SPEI
using Serialization


include("test-ZSI.jl")
nanmaximum(x) = maximum(x[.!isnan.(x)])

wb = deserialize("../data/wb")
x = wb[wb .> 0]

@testset "gamma and loglogistic" begin
  @test fit_gamma(x) == (1.2729799501897179, 17.435758906451866)
  @test fit_logLogistic(x) == (-22.35172530856679, 40.88689577169438, 4.418130441097305)
end

@testset "spei" begin
  @test spei(x).coef == (-22.35172530856679, 40.88689577169438, 4.418130441097305)
  @test spi(x).coef == (1.2729799501897179, 17.435758906451866)
end

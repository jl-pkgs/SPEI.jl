using Statistics, Dates, Test
using Ipaper

@testset "drought_ZSI" begin
  dates = Date(2000):Day(8):Date(2010, 12, 31) |> collect
  A = rand(Float32, 100, 100, length(dates))
  ref = (2001, 2005)
  @time R = drought_ZSI(A, dates; ref)

  by = get_dn.(dates; delta=8)
  grps = unique_sort(by)

  grp = grps[1]
  I = findall(by .== grp)
  I_ref = I[findall(ref[1] .<= year.(dates[I]) .<= ref[2])]

  x = A[1, 1, I]
  x_ref = A[1, 1, I_ref]
  μ = mean(x_ref)
  sd = std(x_ref)
  @test R[1, 1, I] ≈ (x .- μ) ./ sd
end

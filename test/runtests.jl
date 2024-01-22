using Test
using SPEI

nanmaximum(x) = maximum(x[.!isnan.(x)])

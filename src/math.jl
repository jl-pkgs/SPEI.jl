qnorm(p::Real) = quantile(Normal(), p)

pow = ^

# function loggamma(xx::Float64)
#   cof = [76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5]
#   y = x = xx
#   tmp = x + 5.5
#   tmp -= (x + 0.5) * log(tmp)
#   ser = 1.000000000190015
#   for j in 1:6
#     y += 1
#     ser += cof[j] / y
#   end
#   return -tmp + log(2.5066282746310005 * ser / x)
# end

export qnorm

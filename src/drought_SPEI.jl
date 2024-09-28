"""
    drought_SPEI(A::AbstractArray{T,3}, dates;
        ref=(2000, 2020),
        fun_date=get_dn, delta::Int=8,
        progress=true, parallel=true) where {T<:Real}

# Arguments
- `fun`: one of `spei!` and `spi!`
- `nmin`: the minimum number of valid `x_ref`, default 3
  + `spei`: `!isnan(x) = true` as valid
  + `spi`: `x > 0` as valid
"""
function drought_SPEI(A::AbstractArray{T,3}, dates;
  ref=(2000, 2020),
  scale::Int=4, nmin::Int=3,
  fun::Function=spei!,
  fun_date=get_dn, delta::Int=8,
  progress=true, parallel=true) where {T<:Real}

  nlon, nlat, ntime = size(A)
  A2 = similar(A) .* T(NaN)
  accu_scale!(A2, A, scale)

  by = fun_date.(dates; delta)
  grps = unique_sort(by)

  nlon, nlat, ntime = size(A)
  R = zeros(T, nlon, nlat, ntime)

  p = Progress(length(grps))
  @inbounds @par parallel for grp in grps
    progress && next!(p)

    I = findall(by .== grp)
    I_ref = I[findall(ref[1] .<= year.(dates[I]) .<= ref[2])]

    for j = 1:nlon, i = 1:nlat
      _x = @view A2[i, j, I] # 如何处理nan values?
      _x_ref = @view A2[i, j, I_ref]

      z = @view R[i, j, I]
      fun(z, _x, _x_ref; nmin)
    end
  end
  return R
end

function drought_SPEI(x::AbstractVector{T}, dates;
  ref=(2000, 2020),
  scale::Int=4, nmin::Int=3,
  fun::Function=spei!,
  fun_date=get_dn, delta::Int=8) where {T<:Real}

  R = similar(x) .* T(NaN)
  x2 = similar(x) .* T(NaN)
  accu_scale!(x2, x, scale)

  by = fun_date.(dates; delta)
  grps = unique_sort(by)

  @inbounds for grp in grps
    I = findall(by .== grp)
    I_ref = I[findall(ref[1] .<= year.(dates[I]) .<= ref[2])]

    _x = @view x2[I]
    _x_ref = @view x2[I_ref]

    z = @view R[I]
    fun(z, _x, _x_ref; nmin)
  end
  return R
end

function accu_scale!(A2::AbstractArray{T,3}, A::AbstractArray{T,3}, scale::Int) where {T<:Real}
  nlon, nlat, ntime = size(A)
  @inbounds @par for k = scale:ntime
    for j = 1:nlon, i = 1:nlat
      _x = @view A[i, j, k-scale+1:k]
      A2[i, j, k] = nanmean(_x)
    end
  end
  A2[:, :, 1:scale-1] .= T(NaN)
  return A2
end

function accu_scale!(x2::AbstractVector{T}, x::AbstractVector{T}, scale::Int) where {T<:Real}
  ntime = length(x)
  for k = scale:ntime
    _x = @view x[k-scale+1:k]
    x2[k] = nanmean(_x)
  end
  x2[1:scale-1] .= T(NaN)
  return x2
end


export drought_SPEI

"""
    drought_SPEI(A::AbstractArray{T,3}, dates;
        ref=(2000, 2020),
        fun_date=get_dn, delta::Int=8,
        progress=true, parallel=true) where {T<:Real}

Calculate the ZSI (Z-Score Index) of the input data `A` along the time dimension.
> Extremely fast!
"""
function drought_SPEI(A::AbstractArray{T,3}, dates;
  ref=(2000, 2020),
  scale=4,
  fun!::Function=spei!,
  fun_date=get_dn, delta::Int=8,
  progress=true, parallel=true) where {T<:Real}

  nlon, nlat, ntime = size(A)
  A2 = similar(A) .* 0
  @inbounds @par for k = scale:ntime
    for j = 1:nlon, i = 1:nlat
      _x = @view A[i, j, k-scale+1:k]
      A2[i, j, k] = nanmean(_x)
    end
  end
  A2[:, :, 1:scale-1] .= T(NaN)

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
      fun!(z, _x, _x_ref)
    end
  end
  R
end


export drought_SPEI

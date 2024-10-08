"""
    drought_ZSI(A::AbstractArray{T,3}, dates;
        ref=(2000, 2020),
        fun_date=get_dn, delta::Int=8,
        progress=true, parallel=true) where {T<:Real}

Calculate the ZSI (Z-Score Index) of the input data `A` along the time dimension.
> Extremely fast!
"""
function drought_ZSI(A::AbstractArray{T,3}, dates;
  ref=(2000, 2020), nmin::Int=3,
  fun_date=get_dn, delta::Int=8,
  progress=true, parallel=true) where {T<:Real}

  by = fun_date.(dates; delta)
  grps = unique_sort(by)

  nlon, nlat, ntime = size(A)
  R = zeros(T, nlon, nlat, ntime)

  p = Progress(length(grps))
  @inbounds @par parallel for grp in grps
    progress && next!(p)

    I = findall(by .== grp)
    I_ref = I[findall(ref[1] .<= year.(dates[I]) .<= ref[2])]

    for j = 1:nlat, i = 1:nlon
      n_valid = 0
      for k in I_ref
        !isnan(A[i, j, k]) && (n_valid += 1)
      end
      n_valid < nmin && continue

      _x_ref = @view A[i, j, I_ref] # 两次嵌套会导致速度变慢
      μ = nanmean(_x_ref)
      sd = nanstd(_x_ref)

      for _i in I
        R[i, j, _i] = (A[i, j, _i] - μ) / sd
      end
    end
  end
  return R
end


export get_dn, drought_ZSI

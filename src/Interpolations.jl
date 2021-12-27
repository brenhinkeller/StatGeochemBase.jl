## --- Simple linear interpolations

    function _linterp1(x, y, xq::Number)
        knot_index = searchsortedfirst(x, xq, Base.Order.ForwardOrdering()) - 1
        𝔦₋ = min(max(knot_index, firstindex(x)), lastindex(x) - 1)
        𝔦₊ = 𝔦₋ + 1
        x₋, x₊ = x[𝔦₋], x[𝔦₊]
        y₋, y₊ = y[𝔦₋], y[𝔦₊]
        f = (xq - x₋) / (x₊ - x₋)
        return f*y₊ + (1-f)*y₋
    end

    function _linterp1(x, y, xq::AbstractArray; extrapolation=NaN)
        i₁, iₙ = firstindex(x), lastindex(x) - 1
        knot_index = searchsortedfirst_vec(x, xq) .- 1
        T = Base.promote_op(*, eltype(y), Float64)
        yq = similar(y, T, size(xq))
        @turbo for i=1:length(knot_index)
            knot_index[i] = min(max(knot_index[i], i₁), iₙ)
        end
        @turbo for i=1:length(knot_index)
            𝔦₋ = knot_index[i]
            𝔦₊ = 𝔦₋ + 1
            x₋, x₊ = x[𝔦₋], x[𝔦₊]
            y₋, y₊ = y[𝔦₋], y[𝔦₊]
            f = (xq[i] - x₋)/(x₊ - x₋)
            yq[i] = f*y₊ + (1-f)*y₋
        end
        return yq
    end


    # Vectorization-friendly searchsortedfirst implementation from Interpolations.jl
    # https://github.com/JuliaMath/Interpolations.jl
    Base.@propagate_inbounds function searchsortedfirst_exp_left(v, xx, lo, hi)
        for i in 0:4
            ind = lo + i
            ind > hi && return ind
            xx <= v[ind] && return ind
        end
        n = 3
        tn2 = 2^n
        tn2m1 = 2^(n-1)
        ind = lo + tn2
        while ind <= hi
            xx <= v[ind] && return searchsortedfirst(v, xx, lo + tn2 - tn2m1, ind, Base.Order.Forward)
            tn2 *= 2
            tn2m1 *= 2
            ind = lo + tn2
        end
        return searchsortedfirst(v, xx, lo + tn2 - tn2m1, hi, Base.Order.Forward)
    end

    function searchsortedfirst_vec(v::AbstractVector, x::AbstractVector)
        issorted(x) || return searchsortedfirst.(Ref(v), x)
        out = zeros(Int, length(x))
        lo = 1
        hi = length(v)
        @inbounds for i in 1:length(x)
            xx = x[i]
            y = searchsortedfirst_exp_left(v, xx, lo, hi)
            out[i] = y
            lo = min(y, hi)
        end
        return out
    end

## --- Linear interpolation, top-level functions

    # Linear interpolation, finding the interpolated `y` value corresponding to the queried `x` values `xq`
    # Knots (`x`-values) must be sorted
    function linterp1(x::AbstractArray, y::AbstractArray, xq)
        issorted(x) || error("knot-vector `x` must be sorted in increasing order")
        return _linterp1(x, y, xq)
    end
    export linterp1

    # Linear interpolation, finding the interpolated `y` value corresponding to the queried `x` values `xq`
    # Knots will be sorted, along with `y` values, if they are not already
    function linterp1s(x::AbstractArray, y::AbstractArray, xq)
        sI = sortperm(x) # indices to construct sorted array
        return _linterp1(x[sI], y[sI], xq)
    end
    export linterp1s


    # Linearly interpolate vector y at index i, returning outboundsval if outside of bounds
    function linterp_at_index(y::AbstractArray, i::Number, outboundsval=float(eltype(y))(NaN))
        if firstindex(y) <= i < lastindex(y)
            𝔦₋ = floor(Int, i)
            𝔦₊ = 𝔦₋ + 1
            f = i - 𝔦₋
            return f*y[𝔦₊] + (1-f)*y[𝔦₋]
        else
            return outboundsval
        end
    end
    export linterp_at_index



## --- Resize and interpolate arrays of colors

    # Linearly interpolate array of colors at positions xq
    function linterp1(x::AbstractArray, image::AbstractArray{<:Color}, xq)
        # Interpolate red, green, and blue vectors separately
        r_interp = linterp1(x, image .|> c -> c.r, xq)
        g_interp = linterp1(x, image .|> c -> c.g, xq)
        b_interp = linterp1(x, image .|> c -> c.b, xq)
        # Convert back to a color
        return RGB.(r_interp,g_interp,b_interp)
    end

    function resize_colormap(cmap::AbstractArray{<:Color}, n::Integer)
        cNum = length(cmap)
        if n<2
            cmap[1:1]
        else
            linterp1(1:cNum,cmap,collect(range(1,cNum,length=n)))
        end
    end
    export resize_colormap

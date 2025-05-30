## --- Helper functions for fast linear interpolation and extrapolation

    # For parsing kwargs
    staticifsymbol(x::Number) = x
    staticifsymbol(x::StaticSymbol) = x
    function staticifsymbol(x::StaticSymbol{:Linear})
        @warn "extrapolation mode `:Linear` is deprecated, use `:linear` instead"
        static(:linear)
    end
    staticifsymbol(x::Symbol) = staticifsymbol(static(x))

    # Vectorization-friendly searchsortedfirst implementation from Interpolations.jl
    # https://github.com/JuliaMath/Interpolations.jl
    Base.@propagate_inbounds function searchsortedfirst_exp_left(v, xᵢ, lo, hi)
        for i in 0:4
            ind = lo + i
            ind > hi && return ind
            xᵢ <= v[ind] && return ind
        end
        n = 3
        tn2 = 2^n
        tn2m1 = 2^(n-1)
        ind = lo + tn2
        while ind <= hi
            xᵢ <= v[ind] && return searchsortedfirst(v, xᵢ, lo + tn2 - tn2m1, ind, Base.Order.Forward)
            tn2 *= 2
            tn2m1 *= 2
            ind = lo + tn2
        end
        return searchsortedfirst(v, xᵢ, lo + tn2 - tn2m1, hi, Base.Order.Forward)
    end

    function searchsortedfirst_vec!(ix::StridedVector, v::AbstractVector, x::AbstractVector)
        if issorted(x)
            lo = firstindex(v)
            hi = lastindex(v)
            @inbounds for i ∈ eachindex(x, ix)
                y = searchsortedfirst_exp_left(v, x[i], lo, hi)
                ix[i] = y
                lo = min(y, hi)
            end
        else
            ix .= searchsortedfirst.(Ref(v), x)
        end
        return ix
    end

    # Additional methods specializing on specific types
    function searchsortedfirst_vec!(ix::StridedVector, v::AbstractRange, x::AbstractVector)
        lo = firstindex(v)
        hi = lastindex(v) + 1
        δ = first(v) - lo
        λ = step(v)
        @inbounds for i ∈ eachindex(x, ix)
            ix[i] = min(max(ceil(Int, (x[i] - δ)/λ), lo), hi)
        end
        return ix
    end
    function searchsortedfirst_vec!(ix::StridedVector, v::UnitRange, x::AbstractVector)
        lo = firstindex(v)
        hi = lastindex(v) + 1
        δ = first(v) - lo
        @inbounds for i ∈ eachindex(x, ix)
            ix[i] = min(max(ceil(Int, x[i] - δ), lo), hi)
        end
        return ix
    end
    function searchsortedfirst_vec!(ix::StridedVector, v::Base.OneTo, x::AbstractVector)
        lo = firstindex(v)
        hi = lastindex(v) + 1
        @inbounds for i ∈ eachindex(x, ix)
            ix[i] = min(max(ceil(Int, x[i]), lo), hi)
        end
        return ix
    end

## --- 1D linear interpolation

    # Implementation for scalar queries (interpolating a single number)
    Base.@propagate_inbounds _linterp1(x, y, xq::Number, interpolate::T, ::T) where {T<:StaticSymbol} = _linterp1(x, y, xq, interpolate)
    Base.@propagate_inbounds function _linterp1(x, y, xq::Number, interpolate::StaticSymbol, extrapolate::Number)
        if first(x) <= xq <= last(x)
            return _linterp1(x, y, xq, interpolate)
        else
            T = float(eltype(y))
            return T(extrapolate)
        end
    end
    Base.@propagate_inbounds function _linterp1(x, y, xq::Number, interpolate::StaticSymbol, extrapolate::StaticSymbol)
        if first(x) <= xq <= last(x)
            return _linterp1(x, y, xq, interpolate)
        else
            return _linterp1(x, y, xq, extrapolate)
        end
    end
    Base.@propagate_inbounds function _linterp1(x, y, xq::Number, ::StaticSymbol{:linear})
        knot_index = searchsortedfirst(x, xq, Base.Order.ForwardOrdering())
        𝔦₊ = min(max(knot_index, firstindex(x)+1), lastindex(x))
        𝔦₋ = 𝔦₊ - 1
        x₋, x₊ = x[𝔦₋], x[𝔦₊]
        f = (xq - x₋) / (x₊ - x₋)
        return f*y[𝔦₊] + (1-f)*y[𝔦₋]
    end
    Base.@propagate_inbounds function _linterp1(x, y, xq::Number, ::StaticSymbol{:nearest})
        knot_index = searchsortedfirst(x, xq, Base.Order.ForwardOrdering())
        𝔦₊ = min(max(knot_index, firstindex(x)+1), lastindex(x))
        𝔦₋ = 𝔦₊ - 1
        x₋, x₊ = x[𝔦₋], x[𝔦₊]
        f = (xq - x₋) / (x₊ - x₋)
        return if f > 0.5
            y[𝔦₊]
        else
            y[𝔦₋]
        end
    end

    # Forward non-scalar queries to in-place method
    _linterp1(x, y, xq, interpolate::StaticSymbol, extrapolate::Union{Number,StaticSymbol}) = _linterp1!(similar(xq, float(eltype(y))), x, y, xq, interpolate, extrapolate)

    # In-place method
    function _linterp1!(yq, x, y, xq, interpolate::StaticSymbol, extrapolate::Union{Number,StaticSymbol})
        @assert eachindex(xq) == eachindex(yq)
        @inbounds for i in eachindex(xq)
            yq[i] = _linterp1(x, y, xq[i], interpolate, extrapolate)
        end
        return yq
    end

    # In-place methods with optional knot index
    # linear interpolation with constant extrapolation
    function _linterp1k!(yq, knot_index, x, y, xq, ::StaticSymbol{:linear}, extrapolate::Number)
        i₁, iₙ = firstindex(x)+1, lastindex(x)
        searchsortedfirst_vec!(knot_index, x, xq)
        @inbounds for i ∈ eachindex(knot_index)
            𝔦 = knot_index[i]
            if i₁ <= 𝔦 <= iₙ
                𝔦₊ = 𝔦
                𝔦₋ = 𝔦₊ - 1
                x₋, x₊ = x[𝔦₋], x[𝔦₊]
                f = (xq[i] - x₋)/(x₊ - x₋)
                yq[i] = f*y[𝔦₊] + (1-f)*y[𝔦₋]
            elseif first(x) == xq[i]
                yq[i] = first(y)
            else
                yq[i] = extrapolate
            end
        end
        return yq
    end
    # nearest-neighbor interpolation with constant extrapolation
    function _linterp1k!(yq, knot_index, x, y, xq, ::StaticSymbol{:nearest}, extrapolate::Number)
        i₁, iₙ = firstindex(x)+1, lastindex(x)
        searchsortedfirst_vec!(knot_index, x, xq)
        @inbounds for i ∈ eachindex(knot_index)
            𝔦 = knot_index[i]
            if i₁ <= 𝔦 <= iₙ
                𝔦₊ = 𝔦
                𝔦₋ = 𝔦₊ - 1
                x₋, x₊ = x[𝔦₋], x[𝔦₊]
                f = (xq[i] - x₋)/(x₊ - x₋)
                if f > 0.5
                    yq[i] = y[𝔦₊]
                else
                    yq[i] = y[𝔦₋]
                end
            elseif first(x) == xq[i]
                yq[i] = first(y)
            else
                yq[i] = extrapolate
            end
        end
        return yq
    end
    # Linear interpolation and extrapolation
    function _linterp1k!(yq, knot_index, x, y, xq, ::StaticSymbol{:linear}, ::StaticSymbol{:linear})
        i₁, iₙ = firstindex(x)+1, lastindex(x)
        searchsortedfirst_vec!(knot_index, x, xq)
        @inbounds @fastmath for i ∈ eachindex(knot_index)
            𝔦₊ = min(max(knot_index[i], i₁), iₙ)
            𝔦₋ = 𝔦₊ - 1
            x₋, x₊ = x[𝔦₋], x[𝔦₊]
            f = (xq[i] - x₋)/(x₊ - x₋)
            yq[i] = f*y[𝔦₊] + (1-f)*y[𝔦₋]
        end
        return yq
    end
    # Nearest-neighbor interpolation and extrapolation
    function _linterp1k!(yq, knot_index, x, y, xq, ::StaticSymbol{:nearest}, ::StaticSymbol{:nearest})
        i₁, iₙ = firstindex(x)+1, lastindex(x)
        searchsortedfirst_vec!(knot_index, x, xq)
        @inbounds @fastmath for i ∈ eachindex(knot_index)
            𝔦₊ = min(max(knot_index[i], i₁), iₙ)
            𝔦₋ = 𝔦₊ - 1
            x₋, x₊ = x[𝔦₋], x[𝔦₊]
            f = (xq[i] - x₋)/(x₊ - x₋)
            if f > 0.5
                yq[i] = y[𝔦₊]
            else
                yq[i] = y[𝔦₋]
            end
        end
        return yq
    end


    """
    ```julia
    yq = linterp1(x, y, xq; interpolate=static(:linear), extrapolate=interpolate)
    ```
    Simple linear interpolation in one dimension. Given a vector of knots `x`
    and values `y`, find the corresponding `y` values at position(s) `xq`.

    Knots `x` must be sorted in increasing order.

    If the optional keyword argument `extrapolate` is set to `:linear` (default),
    `xq` values outside the range of `x` will be extrapolated using a linear
    extrapolation of the closest two `x`-`y` pairs. Otherwise, if `extrapolate`
    is set to a `Number` (e.g., `0`, or `NaN`), that number will be used
    instead.

    ### Examples
    ```julia
    julia> linterp1(1:10, 1:10, 5.5)
    5.5

    julia> linterp1(1:10, 1:10, 0.5:10.5)
    11-element Vector{Float64}:
      0.5
      1.5
      2.5
      3.5
      4.5
      5.5
      6.5
      7.5
      8.5
      9.5
     10.5
    ```
    """
    function linterp1(x, y, xq; interpolate=static(:linear), extrapolate=interpolate)
        @assert issorted(x) "knot-vector `x` must be sorted in increasing order"
        @assert eachindex(x) == eachindex(y)
        return _linterp1(x, y, xq, staticifsymbol(interpolate), staticifsymbol(extrapolate))
    end
    export linterp1

    """
    ```julia
    linterp1!(yq::AbstractArray, [knot_index], x, y, xq; interpolate=static(:linear), extrapolate=interpolate)
    ```
    In-place variant of `linterp1`.

    An optional temporary working array `knot_index = similar(xq, Int)` may be provided if desired for performance.
    """
    function linterp1!(yq, x, y, xq; interpolate=static(:linear), extrapolate=interpolate)
        @assert issorted(x) "knot-vector `x` must be sorted in increasing order"
        @assert eachindex(x) == eachindex(y)
        @assert eachindex(xq) == eachindex(yq)
        return _linterp1k!(yq, similar(yq, Int), x, y, xq, staticifsymbol(interpolate), staticifsymbol(extrapolate))
    end
    function linterp1!(yq, knot_index::StridedArray{Int}, x, y, xq; interpolate=static(:linear), extrapolate=interpolate)
        @assert issorted(x) "knot-vector `x` must be sorted in increasing order"
        @assert eachindex(x) == eachindex(y)
        @assert eachindex(knot_index) == eachindex(xq) == eachindex(yq)
        return _linterp1k!(yq, knot_index, x, y, xq, staticifsymbol(interpolate), staticifsymbol(extrapolate))
    end
    export linterp1!

    """
    ```julia
    yq = linterp1s(x::AbstractArray, y::AbstractArray, xq; interpolate=static(:linear), extrapolate=interpolate)
    ```
    As as `linterp1` (simple linear interpolation in one dimension), but will sort
    the knots `x` and values `y` pairwise if `x` is not already sorted in increasing order.

    ### Examples
    ```julia
    julia> linterp1s(10:-1:1, 1:10, 5.5)
    5.5

    julia> linterp1s(10:-1:1, 1:10, 0.5:10.5)
    11-element Vector{Float64}:
     10.5
      9.5
      8.5
      7.5
      6.5
      5.5
      4.5
      3.5
      2.5
      1.5
      0.5
    ```
    """
    function linterp1s(x::AbstractArray, y::AbstractArray, xq; interpolate=static(:linear), extrapolate=interpolate)
        @assert eachindex(x) == eachindex(y)
        sI = sortperm(x) # indices to construct sorted array
        return @inbounds _linterp1(x[sI], y[sI], xq, staticifsymbol(interpolate), staticifsymbol(extrapolate))
    end
    export linterp1s

    """
    ```julia
    linterp1s!(yq, [knot_index], x, y, xq; interpolate=static(:linear), extrapolate=interpolate)
    ```
    In-place variant of `linterp1s`.
    Will sort `x` and permute `y` to match, before interpolating at `xq` and storing the result in `yq`.

    An optional temporary working array `knot_index = similar(xq, Int)` may be provided if desired for performance.
    """
    function linterp1s!(yq::StridedArray, x::StridedArray, y::StridedArray, xq; interpolate=static(:linear), extrapolate=interpolate)
        @assert eachindex(x) == eachindex(y)
        @assert eachindex(xq) == eachindex(yq)
        nanargsort!(y, x) # Sort x and permute y to match, in-place
        return _linterp1k!(yq, similar(yq, Int), x, y, xq, staticifsymbol(interpolate), staticifsymbol(extrapolate))
    end
    function linterp1s!(yq::StridedArray, knot_index::StridedArray{Int}, x::StridedArray, y::StridedArray, xq::AbstractArray; interpolate=static(:linear), extrapolate=interpolate)
        @assert eachindex(x) == eachindex(y)
        @assert eachindex(knot_index) == eachindex(yq) == eachindex(xq)
        nanargsort!(y, x) # Sort x and permute y to match, in-place
        return _linterp1k!(yq, knot_index, x, y, xq, staticifsymbol(interpolate), staticifsymbol(extrapolate))
    end
    export linterp1s!

## -- 2D linear interpolation

    function _linterp2(x, y, z::AbstractMatrix, xq, yq, interpolate::StaticSymbol, extrapolate::Union{Number,StaticSymbol})
        # Allocate and fill result
        zq = similar(xq, float(eltype(z)))
        @inbounds for i in eachindex(zq)
            zq[i] = _linterp2(x, y, z, xq[i], yq[i], interpolate, extrapolate)
        end
        return zq
    end
    Base.@propagate_inbounds _linterp2(x, y, z::AbstractMatrix, xq::Number, yq::Number, interpolate::T, ::T) where {T<:StaticSymbol} = _linterp2(x, y, z, xq, yq, interpolate)
    Base.@propagate_inbounds function _linterp2(x, y, z::AbstractMatrix, xq::Number, yq::Number, interpolate::StaticSymbol, extrapolate::Number)
        if first(x) <= xq <= last(x) && first(y) <= yq <= last(y)
            return _linterp2(x, y, z, xq, yq, interpolate)
        else
            T = float(eltype(z))
            return T(extrapolate)
        end
    end
    Base.@propagate_inbounds function _linterp2(x, y, z::AbstractMatrix, xq::Number, yq::Number, interpolate::StaticSymbol, extrapolate::StaticSymbol)
        if first(x) <= xq <= last(x) && first(y) <= yq <= last(y)
            return _linterp2(x, y, z, xq, yq, interpolate)
        else
            return _linterp2(x, y, z, xq, yq, extrapolate)
        end
    end
    Base.@propagate_inbounds function _linterp2(x, y, z::AbstractMatrix, xq::Number, yq::Number, ::StaticSymbol{:nearest})
        x_knot_index = searchsortedfirst(x, xq, Base.Order.ForwardOrdering())
        i₊x = min(max(x_knot_index, firstindex(x)+1), lastindex(x))
        i₋x = i₊x - 1
        x₋, x₊ = x[i₋x], x[i₊x]
        y_knot_index = searchsortedfirst(y, yq, Base.Order.ForwardOrdering()) 
        i₊y = min(max(y_knot_index, firstindex(y)+1), lastindex(y))
        i₋y = i₊y - 1
        y₋, y₊ = y[i₋y], y[i₊y]
        fx = (xq - x₋) / (x₊ - x₋)
        fy = (yq - y₋) / (y₊ - y₋)

        return if fx > 0.5
            if fy > 0.5
                z[i₊y, i₊x]
            else
                z[i₋y, i₊x]
            end
        else
            if fy > 0.5
                z[i₊y, i₋x]
            else
                z[i₋y, i₋x]
            end
        end
    end
    Base.@propagate_inbounds function _linterp2(x, y, z::AbstractMatrix, xq::Number, yq::Number, ::StaticSymbol{:bilinear})
        x_knot_index = searchsortedfirst(x, xq, Base.Order.ForwardOrdering())
        i₊x = min(max(x_knot_index, firstindex(x)+1), lastindex(x))
        i₋x = i₊x - 1
        x₋, x₊ = x[i₋x], x[i₊x]
        y_knot_index = searchsortedfirst(y, yq, Base.Order.ForwardOrdering()) 
        i₊y = min(max(y_knot_index, firstindex(y)+1), lastindex(y))
        i₋y = i₊y - 1
        y₋, y₊ = y[i₋y], y[i₊y]
        z₋₋ = z[i₋y, i₋x]
        z₋₊ = z[i₋y, i₊x]
        z₊₋ = z[i₊y, i₋x]
        z₊₊ = z[i₊y, i₊x]
        fx = (xq - x₋) / (x₊ - x₋)
        fy = (yq - y₋) / (y₊ - y₋)
        z₊ = fx*z₊₊ + (1-fx)*z₊₋
        z₋ = fx*z₋₊ + (1-fx)*z₋₋
        return fy*z₊ + (1-fy)*z₋
    end

    """
    ```julia
    zq = linterp2(x, y, z::AbstractMatrix, xq, yq; interpolate=:bilinear, extrapolate=interpolate)
    ```
    Simple linear interpolation in one dimension. Given vectors of knots `x` and `y`
    and a matrix of values `z`, find the corresponding `z` values at position `xq`,`yq`.

    Knot vectors `x` and `y` must be sorted in increasing order, and must match z
    in dimension, such that `eachindex(x) == axes(z,2)` and `eachindex(y) == axes(z,1)`

    Available methods currently include `:bilinear` (default), and `:nearest` (nearest-neighbor).
    In the `:bilinear` method, `xq`,`yq` pairs will be interpolated (or extrapolated) linearly in `x` 
    and then linearly in `y` (yielding a quadratic result as a whole), based on the 
    closest four z values. In the nearest-neighbor method, the single closest z value will be chosen.
    
    If not otherwise specified, the same method will be used for extrapolation as for interpolation.
    Alternatively, if `extrapolate` is set to a `Number` (e.g., `0`, or `NaN`), that number will be used instead.

    ### Examples
    ```julia
    julia> x = 1:3
     1:3

    julia> y = 1:4
     1:4

    julia> z = y*x'
    4×3 Matrix{Int64}:
     1  2   3
     2  4   6
     3  6   9
     4  8  12

    julia> linterp2(x,y,z,1.5,1.5)
     2.25

    julia> linterp2(x,y,z,2.5,3.5)
     8.75

    julia> linterp2(x,y,z,1,-10)
     -10.0

    julia> linterp2(x,y,z,2,-10)
     -20.0

    julia> linterp2(x,y,z,2,-10, extrapolate=NaN)
     NaN
    ```
    """
    function linterp2(x, y, z::AbstractMatrix, xq, yq; interpolate=static(:bilinear), extrapolate=interpolate)
        @assert issorted(x) "knot-vector `x` must be sorted in increasing order"
        @assert issorted(y) "knot-vector `y` must be sorted in increasing order"
        @assert eachindex(x) == axes(z,2) "Dimensions of `x` must match the horizontal axis of `z`"
        @assert eachindex(y) == axes(z,1) "Dimensions of `y` must match the vertical axis of `z`"
        @assert eachindex(xq) == eachindex(yq) "Dimensions of `xq` and `yq` must match"
        return @inbounds _linterp2(x, y, z, xq, yq, staticifsymbol(interpolate), staticifsymbol(extrapolate))
    end
    export linterp2

## --- Other specialized interpolation methods

    # linearly interpolate vector y at index i, returning outboundsval if outside of bounds
    function linterp_at_index(y::AbstractArray, i::Number, extrapolate=float(eltype(y))(NaN))
        if firstindex(y) <= i < lastindex(y)
            𝔦₋ = floor(Int, i)
            𝔦₊ = 𝔦₋ + 1
            f = i - 𝔦₋
            return f*y[𝔦₊] + (1-f)*y[𝔦₋]
        else
            return extrapolate
        end
    end
    export linterp_at_index


    # linearly interpolate array of colors at positions xq
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

## --- End of File


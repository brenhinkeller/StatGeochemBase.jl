## --- Simple linear interpolations

    function _linterp1(x, y, xq::Number, extrapolate::Symbol)
        @assert extrapolate === :Linear
        knot_index = searchsortedfirst(x, xq, Base.Order.ForwardOrdering()) - 1
        𝔦₋ = min(max(knot_index, firstindex(x)), lastindex(x) - 1)
        𝔦₊ = 𝔦₋ + 1
        x₋, x₊ = x[𝔦₋], x[𝔦₊]
        y₋, y₊ = y[𝔦₋], y[𝔦₊]
        f = (xq - x₋) / (x₊ - x₋)
        return f*y₊ + (1-f)*y₋
    end

    function _linterp1(x, y, xq::Number, extrapolate::Number)
        i₁, iₙ = firstindex(x), lastindex(x) - 1
        knot_index = searchsortedfirst(x, xq, Base.Order.ForwardOrdering()) - 1
        Tₓ = promote_type(eltype(x), eltype(xq))
        T = promote_type(eltype(y), Base.promote_op(/, Tₓ, Tₓ))
        if i₁ <= knot_index <= iₙ
            𝔦₋ = knot_index
            𝔦₊ = 𝔦₋ + 1
            x₋, x₊ = x[𝔦₋], x[𝔦₊]
            y₋, y₊ = y[𝔦₋], y[𝔦₊]
            f = (xq - x₋) / (x₊ - x₋)
            return f*y₊ + (1-f)*y₋
        else
            return T(extrapolate)
        end
    end

    function _linterp1(x, y, xq::AbstractArray, extrapolate)
        Tₓ = promote_type(eltype(x), eltype(xq))
        T = promote_type(eltype(y), Base.promote_op(/, Tₓ, Tₓ))
        yq = similar(xq, T, size(xq))
        _linterp1!(yq, x, y, xq, extrapolate)
    end

    # Allocate knot_index if not provided
    _linterp1!(yq, x, y, xq::AbstractArray, extrapolate) = _linterp1!(yq, ones(Int, length(xq)), x, y, xq::AbstractArray, extrapolate)

    # Linear interpolation with linear extrapolation
    function _linterp1!(yq, knot_index, x::DenseArray, y::DenseArray, xq::AbstractArray, extrapolate::Symbol)
        @assert extrapolate === :Linear
        i₁, iₙ = firstindex(x), lastindex(x) - 1
        searchsortedfirst_vec!(knot_index, x, xq)
        knot_index .-= 1
        @inbounds @fastmath for i ∈ eachindex(knot_index)
            knot_index[i] = min(max(knot_index[i], i₁), iₙ)
        end
        @inbounds @fastmath for i ∈ eachindex(knot_index, xq, yq)
            𝔦₋ = knot_index[i]
            𝔦₊ = 𝔦₋ + 1
            x₋, x₊ = x[𝔦₋], x[𝔦₊]
            y₋, y₊ = y[𝔦₋], y[𝔦₊]
            f = (xq[i] - x₋)/(x₊ - x₋)
            yq[i] = f*y₊ + (1-f)*y₋
        end
        return yq
    end
    # Fallback method
    function _linterp1!(yq, knot_index, x, y, xq::AbstractArray, extrapolate::Symbol)
        @assert extrapolate === :Linear
        i₁, iₙ = firstindex(x), lastindex(x) - 1
        searchsortedfirst_vec!(knot_index, x, xq)
        knot_index .-= 1
        @inbounds for i ∈ eachindex(knot_index)
            knot_index[i] = min(max(knot_index[i], i₁), iₙ)
        end
        @inbounds for i ∈ eachindex(knot_index, xq, yq)
            𝔦₋ = knot_index[i]
            𝔦₊ = 𝔦₋ + 1
            x₋, x₊ = x[𝔦₋], x[𝔦₊]
            y₋, y₊ = y[𝔦₋], y[𝔦₊]
            f = (xq[i] - x₋)/(x₊ - x₋)
            yq[i] = f*y₊ + (1-f)*y₋
        end
        return yq
    end

    # Linear interpolation with constant extrapolation
    function _linterp1!(yq, knot_index, x, y, xq::AbstractArray, extrapolate::Number)
        i₁, iₙ = firstindex(x), lastindex(x) - 1
        searchsortedfirst_vec!(knot_index, x, xq)
        knot_index .-= 1
        @inbounds for i ∈ eachindex(knot_index)
            𝔦 = knot_index[i]
            if i₁ <= 𝔦 <= iₙ
                𝔦₋ = 𝔦
                𝔦₊ = 𝔦₋ + 1
                x₋, x₊ = x[𝔦₋], x[𝔦₊]
                y₋, y₊ = y[𝔦₋], y[𝔦₊]
                f = (xq[i] - x₋)/(x₊ - x₋)
                yq[i] = f*y₊ + (1-f)*y₋
            else
                yq[i] = extrapolate
            end
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

    function searchsortedfirst_vec!(ix::StridedVector, v::AbstractVector, x::AbstractVector)
        @assert firstindex(v) === 1
        if issorted(x)
            lo = 1
            hi = length(v)
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


## --- Linear interpolation, top-level functions


    """
    ```julia
    yq = linterp1(x::AbstractArray, y::AbstractArray, xq; extrapolate=:Linear)
    ```
    Simple linear interpolation in one dimension. Given a vector of knots `x`
    and values `y`, find the corresponding `y` values at position(s) `xq`.

    Knots `x` must be sorted in increasing order.

    If the optional keyword argument `extrapolate` is set to `:Linear` (default),
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
    function linterp1(x::AbstractArray, y::AbstractArray, xq; extrapolate=:Linear)
        issorted(x) || error("knot-vector `x` must be sorted in increasing order")
        return _linterp1(x, y, xq, extrapolate)
    end
    export linterp1

    """
    ```julia
    linterp1!(yq::StridedArray, x::AbstractArray, y::AbstractArray, xq; extrapolate=:Linear, knot_index=ones(Int, length(xq)))
    ```
    In-place variant of `linterp1`.
    """
    function linterp1!(yq::StridedArray, x::AbstractArray, y::AbstractArray, xq; extrapolate=:Linear, knot_index::AbstractVector{Int}=ones(Int, length(xq)))
        issorted(x) || error("knot-vector `x` must be sorted in increasing order")
        return _linterp1!(yq, knot_index, x, y, xq, extrapolate)
    end
    export linterp1!

    """
    ```julia
    yq = linterp1s(x::AbstractArray, y::AbstractArray, xq; extrapolate=:Linear)
    ```
    As as `linterp1` (simple linear interpolation in one dimension), but will sort
    the knots `x` and values `y` pairwise if `x` if not already sorted in
    increasing order.

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
    function linterp1s(x::AbstractArray, y::AbstractArray, xq; extrapolate=:Linear)
        sI = sortperm(x) # indices to construct sorted array
        return _linterp1(x[sI], y[sI], xq, extrapolate)
    end
    export linterp1s

    """
    ```julia
    linterp1s!(yq::StridedArray, x::StridedArray, y::StridedArray, xq; extrapolate=:Linear)
    linterp1s!(yq::StridedArray, knot_index::StridedArray{Int}, x::StridedArray, y::StridedArray, xq::AbstractArray; extrapolate=:Linear)
    ```
    In-place variant of `linterp1s`.
    Will sort `x` and permute `y` to match, before interpolating at `xq` and storing the result in `yq`.

    An optional temporary working array `knot_index = similar(xq, Int)` may be provided to fully eliminate allocations.
    """
    function linterp1s!(yq::StridedArray, x::StridedArray, y::StridedArray, xq; extrapolate=:Linear)
        @assert length(xq) === length(yq)
        @assert eachindex(x) === eachindex(y)
        vsort!(y, x) # Sort x and permute y to match
        return _linterp1!(yq, x, y, xq, extrapolate)
    end
    function linterp1s!(yq::StridedArray, knot_index::StridedArray{Int}, x::StridedArray, y::StridedArray, xq::AbstractArray; extrapolate=:Linear)
        @assert eachindex(knot_index) === eachindex(yq)
        @assert eachindex(x) === eachindex(y)
        @assert length(yq) === length(xq)
        vsort!(y, x) # Sort x and permute y to match
        return _linterp1!(yq, knot_index, x, y, xq, extrapolate)
    end
    export linterp1s!


    # Linearly interpolate vector y at index i, returning outboundsval if outside of bounds
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

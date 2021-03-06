## --- Simple linear interpolations

    function _linterp1(x, y, xq::Number, ::Symbol)
        knot_index = searchsortedfirst(x, xq, Base.Order.ForwardOrdering()) - 1
        ğ¦â = min(max(knot_index, firstindex(x)), lastindex(x) - 1)
        ğ¦â = ğ¦â + 1
        xâ, xâ = x[ğ¦â], x[ğ¦â]
        yâ, yâ = y[ğ¦â], y[ğ¦â]
        f = (xq - xâ) / (xâ - xâ)
        return f*yâ + (1-f)*yâ
    end

    function _linterp1(x, y, xq::AbstractArray, ::Symbol)
        iâ, iâ = firstindex(x), lastindex(x) - 1
        knot_index = searchsortedfirst_vec(x, xq) .- 1
        Tâ = promote_type(eltype(x), eltype(xq))
        T = promote_type(eltype(y), Base.promote_op(/, Tâ, Tâ))
        yq = similar(xq, T, size(xq))
        @inbounds for i â eachindex(knot_index)
            knot_index[i] = min(max(knot_index[i], iâ), iâ)
        end
        @inbounds for i â eachindex(knot_index)
            ğ¦â = knot_index[i]
            ğ¦â = ğ¦â + 1
            xâ, xâ = x[ğ¦â], x[ğ¦â]
            yâ, yâ = y[ğ¦â], y[ğ¦â]
            f = (xq[i] - xâ)/(xâ - xâ)
            yq[i] = f*yâ + (1-f)*yâ
        end
        return yq
    end

    function _linterp1(x, y::AbstractArray{<:AbstractFloat}, xq::AbstractArray, ::Symbol)
        iâ, iâ = firstindex(x), lastindex(x) - 1
        knot_index = searchsortedfirst_vec(x, xq) .- 1
        Tâ = promote_type(eltype(x), eltype(xq))
        T = promote_type(eltype(y), Base.promote_op(/, Tâ, Tâ))
        yq = similar(xq, T, size(xq))
        @turbo for i â eachindex(knot_index)
            knot_index[i] = min(max(knot_index[i], iâ), iâ)
        end
        @turbo for i â eachindex(knot_index)
            ğ¦â = knot_index[i]
            ğ¦â = ğ¦â + 1
            xâ, xâ = x[ğ¦â], x[ğ¦â]
            yâ, yâ = y[ğ¦â], y[ğ¦â]
            f = (xq[i] - xâ)/(xâ - xâ)
            yq[i] = f*yâ + (1-f)*yâ
        end
        return yq
    end

    function _linterp1(x, y, xq::Number, extrapolate::Number)
        iâ, iâ = firstindex(x), lastindex(x) - 1
        knot_index = searchsortedfirst(x, xq, Base.Order.ForwardOrdering()) - 1
        Tâ = promote_type(eltype(x), eltype(xq))
        T = promote_type(eltype(y), Base.promote_op(/, Tâ, Tâ))
        if iâ <= knot_index <= iâ
            ğ¦â = knot_index
            ğ¦â = ğ¦â + 1
            xâ, xâ = x[ğ¦â], x[ğ¦â]
            yâ, yâ = y[ğ¦â], y[ğ¦â]
            f = (xq - xâ) / (xâ - xâ)
            return f*yâ + (1-f)*yâ
        else
            return T(extrapolate)
        end
    end

    function _linterp1(x, y, xq::AbstractArray, extrapolate::Number)
        iâ, iâ = firstindex(x), lastindex(x) - 1
        knot_index = searchsortedfirst_vec(x, xq) .- 1
        Tâ = promote_type(eltype(x), eltype(xq))
        T = promote_type(eltype(y), Base.promote_op(/, Tâ, Tâ))
        yq = similar(xq, T, size(xq))
        @inbounds for i â eachindex(knot_index)
            ğ¦ = knot_index[i]
            if iâ <= ğ¦ <= iâ
                ğ¦â = ğ¦
                ğ¦â = ğ¦â + 1
                xâ, xâ = x[ğ¦â], x[ğ¦â]
                yâ, yâ = y[ğ¦â], y[ğ¦â]
                f = (xq[i] - xâ)/(xâ - xâ)
                yq[i] = f*yâ + (1-f)*yâ
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


    # Linearly interpolate vector y at index i, returning outboundsval if outside of bounds
    function linterp_at_index(y::AbstractArray, i::Number, extrapolate=float(eltype(y))(NaN))
        if firstindex(y) <= i < lastindex(y)
            ğ¦â = floor(Int, i)
            ğ¦â = ğ¦â + 1
            f = i - ğ¦â
            return f*y[ğ¦â] + (1-f)*y[ğ¦â]
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

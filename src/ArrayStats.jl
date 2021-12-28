## --- Sorting and counting array elements

    """
    ```julia
    n = count_unique!(A)
    ```
    Sort the array `A` in-place (if not already sorted), move unique elements to
    the front, and return the number of unique elements found.

    `A[1:count_unique!(A)]` should return an array equivalent to `unique(A)`
    """
    function count_unique!(A)
        issorted(A) || sort!(A)
        n = 1
        last = A[1]
        @inbounds for i=2:length(A)
            if A[i] != last
                n += 1
                last = A[n] = A[i]
            end
        end
        return n
    end
    export count_unique!


## --- Convert between bin centers and bin edges

    """
    ```julia
    cntr(edges::AbstractArray)
    ```
    Given an array of bin edges, return a corresponding vector of bin centers
    """
    function cntr(edges::AbstractArray)
        centers = (edges[1:end-1] + edges[2:end]) ./ 2
        return centers
    end
    export cntr


## --- Searching arrays

    """
    ```julia
    findmatches(source, target)
    ```
    Return the index of the first value in `target` (if any) that is equal to
    a given value in `source` for each value in `source`; else 0.
    """
    function findmatches(source, target)
        index = Array{Int64}(undef, size(source))
        return findmatches!(index, source, target)
    end
    function findmatches!(index, source, target)
        # Loop through source and find first match for each (if any)
        @inbounds for i = 1:length(index)
            for j = 1:length(target)
                index[i] = 0
                if isequal(source[i], target[j])
                    index[i] = j
                    break
                end
            end
        end
        return index
    end
    export findmatches, findmatches!

    """
    ```julia
    findclosest(source, target)
    ```
    Return the index of the numerically closest value in the indexable collection
    `target` for each value in `source`.
    If muliple values are equally close, the first one is used
    """
    function findclosest(source, target)
        index = Array{Int64}(undef, size(source))
        return findclosest!(index, source, target)
    end
    function findclosest!(index, source, target)
        # Find closest (numerical) match in target for each value in source
        @inbounds for i = 1:length(source)
            d = abs(target[1] - source[i])
            index[i] = 1
            for j = 2:length(target)
                d_prop = abs(target[j] - source[i])
                if d_prop < d
                    d = d_prop
                    index[i] = j
                end
            end
        end
        return index
    end
    export findclosest, findclosest!

    """
    ```julia
    findclosestbelow(source, target)
    ```
    Return the index of the nearest value of the indexable collection `target`
    that is less than (i.e., "below") each value in `source`.
    If no such target values exist in `target`, returns an index of 0.
    """
    function findclosestbelow(source, target)
        index = Array{Int64}(undef, size(source))
        return findclosestbelow!(index, source, target)
    end
    function findclosestbelow!(index, source, target)
        @inbounds for i = 1:length(source)
            index[i] = d = j = 0
            while j < length(target)
                j += 1
                if target[j] < source[i]
                    d = source[i] - target[j]
                    index[i] = j
                    break
                end
            end
            while j < length(target)
                j += 1
                if target[j] < source[i]
                    d_prop = source[i] - target[j]
                    if d_prop < d
                        d = d_prop
                        index[i] = j
                    end
                end
            end
        end
        return index
    end
    export findclosestbelow, findclosestbelow!

    """
    ```julia
    findclosestabove(source, target)
    ```
    Return the index of the nearest value of the indexable collection `target`
    that is greater than (i.e., "above") each value in `source`.
    If no such values exist in `target`, returns an index of 0.
    """
    function findclosestabove(source, target)
        index = Array{Int64}(undef, size(source))
        return findclosestabove!(index,source,target)
    end
    function findclosestabove!(index, source, target)
        @inbounds for i = 1:length(source)
            index[i] = d = j = 0
            while j < length(target)
                j += 1
                if target[j] > source[i]
                    d = target[j] - source[i]
                    index[i] = j
                    break
                end
            end
            while j < length(target)
                j += 1
                if target[j] > source[i]
                    d_prop = target[j] - source[i]
                    if d_prop < d
                        d = d_prop
                        index[i] = j
                    end
                end
            end
        end
        return index
    end
    export findclosestabove, findclosestabove!

    """
    ```julia
    findnth(t::AbstractArray{Bool}, n::Integer)
    ```
    Return the index of the `n`th true value of `t`, else length(`t`)
    """
    function findnth(t::AbstractArray{Bool}, n::Integer)
        N = 0
        @inbounds for i ∈ eachindex(t)
            if t[i]
                N += 1
            end
            if N == n
                return i
            end
        end
        return length(t)
    end
    export findnth

## --- String matching

    # if ~ @isdefined contains
    #     """
    #     ```julia
    #     contains(haystack, needle)
    #     ```
    #
    #     Converts both `haystack` and `needle` to strings (if not already strings)
    #     and checks whether `string(haystack)` contains `string(needle)`.
    #     """
    #     contains(haystack::AbstractString, needle::Union{AbstractString,Regex,AbstractChar}) = occursin(needle, haystack)
    #     contains(haystack, needle) = occursin(string(needle), string(haystack))
    #     export contains
    # end

    """
    ```julia
    containsi(haystack, needle)
    ```

    Converts both `haystack` and `needle` to strings and checks whether
    `string(haystack)` contains `string(needle)`, ignoring case.
    """
    containsi(haystack::AbstractString, needle::Union{AbstractString,AbstractChar}) = occursin(lowercase(needle), lowercase(haystack))
    containsi(haystack, needle) = occursin(lowercase(string(needle)), lowercase(string(haystack)))
    export containsi


## --- Drawing a pseudorandom array from a numerically specified distribution

    """
    ```julia
    x = draw_from_distribution(dist::AbstractArray{<:AbstractFloat}, n::Integer)
    ```
    Draw `n` random floating point numbers from a continuous probability distribution
    specified by a vector `dist` defining the PDF curve thereof.
    """
    function draw_from_distribution(dist::AbstractArray{<:AbstractFloat}, n::Integer)
        x = Array{eltype(dist)}(undef, n)
        draw_from_distribution!(x, dist)
        return x
    end
    export draw_from_distribution

    """
    ```julia
    draw_from_distribution!(dist::AbstractArray{<:AbstractFloat}, x::Array{<:AbstractFloat})
    ```
    Fill an existing variable `x` with random floating point numbers drawn from
    a continuous probability distribution specified by a vector `dist`
    defining the PDF curve thereof.
    """
    function draw_from_distribution!(x::Array{<:AbstractFloat}, dist::AbstractArray{<:AbstractFloat})
        # Fill the array x with random numbers from the distribution 'dist'
        dist_ymax = maximum(dist)
        dist_xmax = prevfloat(length(dist) - 1.0)

        @inbounds for i ∈ eachindex(x)
            while true
                # Pick random x value
                rx = rand(eltype(x)) * dist_xmax
                # Interpolate corresponding distribution value
                f = floor(Int,rx)
                y = dist[f+2]*(rx-f) + dist[f+1]*(1-(rx-f))
                # See if x value is accepted
                ry = rand(Float64) * dist_ymax
                if (y > ry)
                    x[i] = rx / dist_xmax
                    break
                end
            end
        end
    end
    export draw_from_distribution!


## --- Numerically integrate a 1-d distribution

    """
    ```julia
    trapezoidalquadrature(edges, values)
    ```
    Add up the area under a curve with y positions specified by a vector of `values`
    and x positions specfied by a vector of `edges` using trapezoidal integration.
    Bins need not be evenly spaced, though it helps (integration will be faster
    if `edges` are specified as an AbstractRange).
    """
    function trapezoidalquadrature(edges::AbstractRange, values::AbstractArray)
        result = zero(eltype(values))
        @turbo for i = (firstindex(edges)+1):lastindex(edges)
            result += values[i-1]+values[i]
        end
        dx = (edges[end]-edges[1])/(length(edges) - 1)
        return result * dx / 2
    end
    function trapezoidalquadrature(edges::AbstractArray, values::AbstractArray)
        result = zero(promote_type(eltype(edges), eltype(values)))
        @turbo for i = (firstindex(edges)+1):lastindex(edges)
            result += (values[i-1] + values[i]) * (edges[i] - edges[i-1])
        end
        return result / 2
    end
    export trapezoidalquadrature
    trapz = trapezoidalquadrature
    export trapz

    """
    ```julia
    midpointquadrature(bincenters, values)
    ```
    Add up the area under a curve with y positions specified by a vector of `values`
    and x positions specfied by a vector of `bincenters` using midpoint integration.
    """
    function midpointquadrature(bincenters::AbstractRange, values::AbstractArray)
        sum(values) * (last(bincenters)-first(bincenters)) / (length(bincenters) - 1)
    end
    export midpointquadrature


## --- End of File

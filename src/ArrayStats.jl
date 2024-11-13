## --- To make arrays with messy types better behaved

    """
    ```julia
    unionize(x::AbstractVector)
    ```
    Turn an array with possibly abstract element type into one with 
    `eltype` equal to a Union of all types of elements in the array.

    Always returns a copy, even if `x` is already unionized.

    ### Examples
    ```julia
    julia> a = Any[false, 0, 1.0]
    3-element Vector{Any}:
     false
         0
         1.0
    
    julia> unionize(a)
    3-element Vector{Union{Bool, Float64, Int64}}:
     false
         0
         1.0
    ```
    """
    function unionize(x::AbstractVector)
        types = unique(typeof.(x))
        if length(types) > 1
            unionized = similar(x, Union{types...})
        else
            unionized = similar(x, only(types))
        end
        unionized .= x
    end
    unionize(x::AbstractRange) = copy(x) # Exemption for ranges, which should probably alway have concrete eltype already
    export unionize


## --- To avoid allocations when indexing by a vector of Booleans


    """
    ```julia
    copyat!(dest, src, tₛ::AbstractVector{Bool})
    ```
    Copy from src to dest when tₛ is true. Equivalent to `dest .= src[tₛ]`, but without inducing allocations.

    See also `reversecopyat!`
    """
    function copyat!(dest::DenseArray, src, tₛ::AbstractVector{Bool})
        @assert eachindex(src) == eachindex(tₛ)
        iₙ = firstindex(dest)
        iₗ = lastindex(dest)
        @inbounds for iₛ in eachindex(src)
            if tₛ[iₛ]
                dest[iₙ] = src[iₛ]
                iₙ += 1
                iₙ > iₗ && break
            end
        end
        return dest
    end
    export copyat!

    """
    ```julia
    reversecopyat!(dest, src, tₛ::AbstractVector{Bool})
    ```
    As `copyat!`, but also reverse the order of stored elements. 
    
    Equivalent to `dest .= reverse(src[tₛ])`, but without inducing allocations.
    """
    function reversecopyat!(dest::DenseArray, src, tₛ::AbstractVector{Bool})
        @assert eachindex(src) == eachindex(tₛ)
        i₀ = firstindex(dest)
        iₙ = lastindex(dest)
        @inbounds for iₛ in eachindex(src)
            if tₛ[iₛ]
                dest[iₙ] = src[iₛ]
                iₙ -= 1
                iₙ < i₀ && break
            end
        end
        return dest
    end
    export reversecopyat!

## --- Sorting and counting array elements

    """
    ```julia
    n = count_unique!(A)
    ```
    Sort the array `A` in-place (if not already sorted), move unique elements to
    the front, and return the number of unique elements found.

    `A[1:count_unique!(A)]` should return an array equivalent to `unique(A)`

    ### Examples
    ```julia
    julia> A = rand(1:5, 10)
    10-element Vector{Int64}:
     4
     4
     2
     3
     3
     4
     1
     5
     1
     2

    julia> A = rand(1:5, 7)
    7-element Vector{Int64}:
     1
     1
     4
     3
     1
     1
     4

    julia> n = count_unique!(A)
    3

    julia> A
    7-element Vector{Int64}:
     1
     3
     4
     1
     3
     4
     4

    julia> A[1:n]
    3-element Vector{Int64}:
     1
     3
     4
    ```
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
    cntr(edges::Collection)
    ```
    Given an array of bin edges, return a corresponding vector of bin centers

    ### Examples
    ```julia
    julia> cntr(1:10)
     1.5:1.0:9.5
    ```
    """
    function cntr(edges::Collection)
        centers = (edges[1:end-1] + edges[2:end]) ./ 2
        return centers
    end
    export cntr

## -- Other utilities for working with binned data and histograms

    """
    ```julia
    stepify(x::AbstractVector)
    ```
    Duplicate every element of an array.
    
    Together with `stepifyedges`, allows for convenient plotting of 
    a histogram of bin counts and bin edges as a stepped line.
    ### Examples
    ```julia
    julia> stepify([1,3,2])
    6-element Vector{Int64}:
     1
     1
     3
     3
     2
     2

    julia> plot(stepifyedges(binedges), stepify(bincounts))
    # Returns a stepped line plot of a histogram
    ```
    """
    stepify(x::AbstractVector) = vec(vcat(x', x'))
    export stepify

    """
    ```julia
    stepifyedges(x::AbstractVector)
    ```
    Duplicate every element of an array except for the first and last.
    
    Together with `stepify`, allows for convenient plotting of 
    a histogram of bin counts and bin edges as a stepped line.
    ### Examples
    ```julia
    julia> stepifyedges([1,2,3,4])
    6-element Vector{Int64}:
     1
     2
     2
     3
     3
     4

    julia> plot(stepifyedges(binedges), stepify(bincounts))
    # Returns a stepped line plot of a histogram
    ```
    """
    stepifyedges(x::AbstractVector) = vec(vcat(x[1:end-1]', x[2:end]'))
    export stepifyedges

## --- Searching arrays

    """
    ```julia
    findmatches(needles, haystack)
    ```
    Return the linear index of the first value in `haystack` (if any) that is equal
    to a given value in `needles` for each value in `needles`; else 0.

    ### Examples
    ```julia
    julia> findmatches([3,5],1:10)
    2-element Vector{Int64}:
     3
     5
    ```
    """
    function findmatches(needle, haystack::Collection)
        @inbounds for j ∈ eachindex(haystack)
            if isequal(needle, haystack[j])
                return j
            end
        end
        return 0
    end
    findmatches(needles::Collection, haystack::Collection) = findmatches!(fill(0, length(needles)), needles, haystack)
    function findmatches!(index::DenseArray, needles, haystack)
        @assert eachindex(index) == eachindex(needles)
        # Loop through needles and find first match for each (if any)
        @inbounds for i ∈ eachindex(index)
            index[i] = 0
            for j ∈ eachindex(haystack)
                if isequal(needles[i], haystack[j])
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
    findclosest(needles, haystack)
    ```
    Return the index of the numerically closest value in the indexable collection
    `haystack` for each value in `needles`.
    If muliple values are equally close, the first one is used

    ### Examples
    ```julia
    julia> findclosest(3.4, 1:10)
    3

    julia> findclosest(3:4, 1:10)
    2-element Vector{Int64}:
     3
     4
    ```
    """
    function findclosest(needle, haystack::Collection)
        if issorted(haystack)
            𝔦ₛ = searchsortedfirst(haystack, needle)
            𝔦₊ = min(𝔦ₛ, lastindex(haystack))
            𝔦₋ = max(𝔦ₛ-1, firstindex(haystack))
            index = if 𝔦₊ != 𝔦₋ && abs(haystack[𝔦₊]-needle) > abs(haystack[𝔦₋]-needle) 
                𝔦₋
            else
                𝔦₊
            end
        elseif issorted(haystack, rev=true)
            𝔦ₛ = searchsortedfirst(haystack, needle, rev=true)
            𝔦₊ = min(𝔦ₛ, lastindex(haystack))
            𝔦₋ = max(𝔦ₛ-1, firstindex(haystack))
            index = if 𝔦₊ != 𝔦₋ && abs(haystack[𝔦₊]-needle) > abs(haystack[𝔦₋]-needle) 
                𝔦₋
            else
                𝔦₊
            end
        else
            δ = abs(first(haystack) - needle)
            index = firstindex(haystack)
            @inbounds for j ∈ Iterators.drop(eachindex(haystack),1)
                δₚ = abs(haystack[j] - needle)
                if δₚ < δ
                    δ = δₚ
                    index = j
                end
            end
        end
        return index
    end
    findclosest(needles::Collection, haystack::Collection) = findclosest!(fill(0, length(needles)), needles, haystack)
    function findclosest!(index::DenseArray, needles, haystack)
        @assert eachindex(index) == eachindex(needles)
        # Find closest (numerical) match in haystack for each value in needles
        if issorted(haystack)
            @inbounds for i ∈ eachindex(needles)
                𝔦ₛ = searchsortedfirst(haystack, needles[i])
                𝔦₊ = min(𝔦ₛ, lastindex(haystack))
                𝔦₋ = max(𝔦ₛ-1, firstindex(haystack))
                if 𝔦₊ != 𝔦₋ && abs(haystack[𝔦₊]-needles[i]) > abs(haystack[𝔦₋]-needles[i]) 
                    index[i] = 𝔦₋
                else
                    index[i] = 𝔦₊
                end
            end
        elseif issorted(haystack, rev=true)
            @inbounds for i ∈ eachindex(needles)
                𝔦ₛ = searchsortedfirst(haystack, needles[i], rev=true)
                𝔦₊ = min(𝔦ₛ, lastindex(haystack))
                𝔦₋ = max(𝔦ₛ-1, firstindex(haystack))
                if 𝔦₊ != 𝔦₋ && abs(haystack[𝔦₊]-needles[i]) > abs(haystack[𝔦₋]-needles[i]) 
                    index[i] = 𝔦₋
                else
                    index[i] = 𝔦₊
                end
            end
        else
            @inbounds for i ∈ eachindex(needles)
                δ = abs(first(haystack) - needles[i])
                index[i] = firstindex(haystack)
                for j ∈ Iterators.drop(eachindex(haystack),1)
                    δₚ = abs(haystack[j] - needles[i])
                    if δₚ < δ
                        δ = δₚ
                        index[i] = j
                    end
                end
            end

        end
        return index
    end
    export findclosest, findclosest!


    """
    ```julia
    findclosestbelow(needles, haystack)
    ```
    Return the index of the nearest value of the indexable collection `haystack`
    that is less than (i.e., "below") each value in `needles`.
    If no such haystack values exist, returns `firstindex(haystack)-1`.

    ### Examples
    ```julia
    julia> findclosestbelow(3.5, 1:10)
    3
    
    julia> findclosestbelow(3:4, 1:10)
    2-element Vector{Int64}:
     2
     3    
    ```
    """
    function findclosestbelow(needle, haystack::Collection)
        if issorted(haystack)
            index = searchsortedfirst(haystack, needle) - 1
        elseif issorted(haystack, rev=true)
            index = searchsortedlast(haystack, needle, rev=true) + 1
            index > lastindex(haystack) && (index = firstindex(haystack)-1)
        else
            δ = needle - first(haystack)
            index = j = firstindex(haystack) - 1
            while j < lastindex(haystack)
                j += 1
                if haystack[j] < needle
                    δ = needle - haystack[j]
                    index = j
                    break
                end
            end
            while j < lastindex(haystack)
                j += 1
                if haystack[j] < needle
                    δₚ = needle - haystack[j]
                    if δₚ < δ
                        δ = δₚ
                        index = j
                    end
                end
            end
        end
        return index
    end
    findclosestbelow(needles::Collection, haystack::Collection) = findclosestbelow!(fill(0, length(needles)), needles, haystack)
    function findclosestbelow!(index::DenseArray, needles, haystack)
        if issorted(haystack)
            @inbounds for i ∈ eachindex(needles)
                index[i] = searchsortedfirst(haystack, needles[i]) - 1
            end
        elseif issorted(haystack, rev=true)
            @inbounds for i ∈ eachindex(needles)
                index[i] = searchsortedlast(haystack, needles[i], rev=true) + 1
                index[i] > lastindex(haystack) && (index[i] = firstindex(haystack)-1)
            end
        else
            ∅ = firstindex(haystack) - 1
            δ = first(needles) - first(haystack)
            @inbounds for i ∈ eachindex(needles)
                index[i] = j = ∅
                while j < lastindex(haystack)
                    j += 1
                    if haystack[j] < needles[i]
                        δ = needles[i] - haystack[j]
                        index[i] = j
                        break
                    end
                end
                while j < lastindex(haystack)
                    j += 1
                    if haystack[j] < needles[i]
                        δₚ = needles[i] - haystack[j]
                        if δₚ < δ
                            δ = δₚ
                            index[i] = j
                        end
                    end
                end
            end
        end
        return index
    end
    export findclosestbelow, findclosestbelow!


    """
    ```julia
    findclosestabove(needles, haystack)
    ```
    Return the index of the nearest value of the indexable collection `haystack`
    that is greater than (i.e., "above") each value in `needles`.
    If no such values exist, returns `lastindex(haystack)+1`.

    ### Examples
    ```julia
    julia> findclosestabove(3.5, 1:10)
    4

    julia> findclosestabove(3:4, 1:10)
    2-element Vector{Int64}:
     4
     5
    ```
    """
    function findclosestabove(needle, haystack::Collection)
        if issorted(haystack)
            index = searchsortedlast(haystack, needle) + 1
        elseif issorted(haystack, rev=true)
            index = searchsortedfirst(haystack, needle, rev=true) - 1
            index < firstindex(haystack) && (index = lastindex(haystack)+1)
        else
            δ = needle - first(haystack)
            index = j = lastindex(haystack) + 1
            while j > firstindex(haystack)
                j -= 1
                if haystack[j] > needle
                    δ = haystack[j] - needle
                    index = j
                    break
                end
            end
            while j > firstindex(haystack)
                j -= 1
                if haystack[j] > needle
                    δₚ = haystack[j] - needle
                    if δₚ < δ
                        δ = δₚ
                        index = j
                    end
                end
            end
        end
        return index
    end
    findclosestabove(needles::Collection, haystack::Collection) = findclosestabove!(fill(0, length(needles)), needles, haystack)
    function findclosestabove!(index::DenseArray, needles, haystack)
        if issorted(haystack)
            @inbounds for i ∈ eachindex(needles)
                index[i] = searchsortedlast(haystack, needles[i]) + 1
            end
        elseif issorted(haystack, rev=true)
            @inbounds for i ∈ eachindex(needles)
                index[i] = searchsortedfirst(haystack, needles[i], rev=true) - 1
                index[i] < firstindex(haystack) && (index[i] = lastindex(haystack)+1)
            end
        else
            ∅ = lastindex(haystack) + 1
            δ = first(needles) - first(haystack)
            @inbounds for i ∈ eachindex(needles)
                index[i] = j = ∅
                while j > firstindex(haystack)
                    j -= 1
                    if haystack[j] > needles[i]
                        δ = haystack[j] - needles[i]
                        index[i] = j
                        break
                    end
                end
                while j > firstindex(haystack)
                    j -= 1
                    if haystack[j] > needles[i]
                        δₚ = haystack[j] - needles[i]
                        if δₚ < δ
                            δ = δₚ
                            index[i] = j
                        end
                    end
                end
            end
        end
        return index
    end
    export findclosestabove, findclosestabove!


    """
    ```julia
    findnth(t::Collection{Bool}, n::Integer)
    ```
    Return the index of the `n`th true value of `t`, else length(`t`)

    ### Examples
    ```julia
    julia> t = rand(Bool,5)
    5-element Vector{Bool}:
     1
     1
     0
     1
     1

    julia> findnth(t, 3)
    4
    ```
    """
    function findnth(t::Collection{Bool}, n::Integer)
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


    """
    ```julia
    findclosestunequal(x::Collection, i::Integer)
    ```
    Return the index of the closest index `n` to `i` for which `x[n] != x[i]`,
    or `i` if no unequal values of `x` are found.

    ### Examples
    ```julia
    julia> x = [1, 2, 2, 3, 4]
    5-element Vector{Int64}:
     1
     2
     2
     3
     4
    
    julia> findclosestunequal(x, 2)
    1
    
    julia> findclosestunequal(x, 3)
    4    
    ```
    """
    function findclosestunequal(x::Collection, i::Int)
        xᵢ = x[i]
        for offset = 1:(length(x)-1)
            l = i - offset
            if l >= firstindex(x)
                (x[l] == xᵢ) || return l
            end
            u = i + offset
            if u <= lastindex(x)
                (x[u] == xᵢ) || return u
            end
        end
        return i
    end
    export findclosestunequal

## --- String matching


    """
    ```julia
    containsi(haystack, needle)
    ```
    Converts both `haystack` and `needle` to strings and checks whether
    `string(haystack)` contains `string(needle)`, ignoring case.

    ### Examples
    ```julia
    julia> containsi("QuickBrownFox", "brown")
    true
    ```
    """
    containsi(haystack::AbstractString, needle::Union{AbstractString,AbstractChar}) = occursin(lowercase(needle), lowercase(haystack))
    containsi(haystack, needle) = occursin(lowercase(string(needle)), lowercase(string(haystack)))
    export containsi


## --- Drawing a pseudorandom array from a numerically specified distribution

    """
    ```julia
    draw_from_distribution(dist::Collection{AbstractFloat}, n::Integer)
    ```
    Draw `n` random floating point numbers from a continuous probability distribution
    specified by a collection `dist` defining the PDF curve thereof.

    ### Examples
    ```julia
    julia> draw_from_distribution([0,1,2,1,0.], 7)
    7-element Vector{Float64}:
     0.5271744125470383
     0.6624591724796276
     0.7737643383545575
     0.9603780726501608
     0.7772477857811155
     0.8307248435614027
     0.6351766227803024    
    ```
    """
    function draw_from_distribution(dist::Collection{AbstractFloat}, n::Integer)
        x = Array{eltype(dist)}(undef, n)
        draw_from_distribution!(x, dist)
        return x
    end
    export draw_from_distribution

    """
    ```julia
    draw_from_distribution!(x::DenseArray{<:AbstractFloat}, dist::Collection{AbstractFloat})
    ```
    Fill an existing variable `x` with random floating point numbers drawn from
    a continuous probability distribution specified by a vector `dist`
    defining the PDF curve thereof.
    """
    function draw_from_distribution!(x::DenseArray{<:AbstractFloat}, dist::Collection{AbstractFloat})
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

    ### Examples
    ```julia
    julia> trapezoidalquadrature(0:0.1:10, 0:0.1:10)
    50.0
    ```
    """
    function trapezoidalquadrature(edges::AbstractRange, values::Collection)
        @assert eachindex(edges)==eachindex(values)
        result = zero(eltype(values))
        @inbounds @fastmath for i ∈ (firstindex(edges)+1):lastindex(edges)
            result += values[i-1]+values[i]
        end
        dx = (edges[end]-edges[1])/(length(edges) - 1)
        return result * dx / 2
    end
    function trapezoidalquadrature(edges::Collection, values::Collection)
        @assert eachindex(edges)==eachindex(values)
        result = zero(promote_type(eltype(edges), eltype(values)))
        @inbounds @fastmath for i ∈ (firstindex(edges)+1):lastindex(edges)
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

    ### Examples
    ```julia
    julia> midpointquadrature(0:0.1:10, 0:0.1:10)
    50.5
    ```
    """
    function midpointquadrature(bincenters::AbstractRange, values::Collection)
        @assert eachindex(bincenters)==eachindex(values)
        sum(values) * (last(bincenters)-first(bincenters)) / (length(bincenters) - 1)
    end
    export midpointquadrature


## --- End of File

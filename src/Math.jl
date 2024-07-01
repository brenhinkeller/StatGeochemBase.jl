## --- Dealing with different number representations

    """
    ```julia
    nearest(T, x)
    ```
    Convert `x` to the nearest representable value in type T, rounding if inexact

    ### Examples
    ```julia
    julia> nearest(Int, 1234.56)
    1235

    julia> nearest(Int, Inf)
    9223372036854775807

    julia> nearest(Int, -Inf)
    -9223372036854775808
    ````
    """
    function nearest(::Type{T}, x) where T <: Number
        if x > typemax(T)
            typemax(T)
        elseif x < typemin(T)
            typemin(T)
        else
            T(x)
        end
    end
    function nearest(::Type{T}, x) where T <: Integer
        if x > typemax(T)
            typemax(T)
        elseif x < typemin(T)
            typemin(T)
        else
            round(T, x)
        end
    end
    export nearest

## --- Determining reported precision of numbers

    # Convert size to decimal precision
    maxdigits(T::Type) = ceil(Int, sizeof(T)*2.408239965311849)
    # Special cases
    maxdigits(::Type{BigFloat}) = 78
    maxdigits(::Type{Float64}) = 16
    maxdigits(::Type{Float32}) = 8
    maxdigits(::Type{Float16}) = 4
    maxdigits(::Type{Int64}) = 19
    maxdigits(::Type{Int32}) = 10
    maxdigits(::Type{Int16}) = 5
    maxdigits(::Type{Int8}) = 3
    maxdigits(::Type{Bool}) = 1

    """
    ```julia
    sigdigits(d)
    ```
    Determine the number of decimal significant figures of a number `d`.

    ### Examples
    ```julia
    julia> sigdigits(1000)
    1
    
    julia> sigdigits(1001)
    4
    
    julia> sigdigits(1001.1245)
    8
    ```
    """
    function sigdigits(d::T) where T <: Number
        n = 0
        isfinite(d) || return n
        rtol = 10.0^-maxdigits(T)
        while n < maxdigits(T)
            isapprox(d, round(d, sigdigits=n); rtol) && return n
            n += 1
        end
        return n 
    end
    sigdigits(d::Irrational) = Inf
    export sigdigits


    """
    ```julia
    leastsigfig(d)
    ```
    Return the order of magnitude of the least significant decimal digit of a number `d`.

    ### Examples
    ```julia
    julia> leastsigfig(1000)
    1000.0
    
    julia> leastsigfig(1001)
    1.0
    
    julia> leastsigfig(1001.1234)
    0.0001
    ```
    """
    function leastsigfig(d)
        iszero(d) && return 1.0*d
        isfinite(d) || return 1.0*d
        10.0^(floor(Int, log10(abs(d)))-sigdigits(d)+1)
    end
    export leastsigfig



## --- Fast inverse square-root

    """
    ```julia
    fast_inv_sqrt(x)
    ```
    The infamous fast inverse square root of `x`, in 32 and 64 bit versions.
    Can be up to 10x faster than base `1/sqrt(x)`, though with nontrivial loss
    of precision. The implementations here are good to about 4 ppm.
    """
    function fast_inv_sqrt(x::Float64)
        x_2 = 0.5 * x
        result = Base.sub_int(9.603007803048109e153, Base.lshr_int(x,1)) # Floating point magic
        result *= ( 1.5 - (x_2 * result * result )) # Newton's method
        result *= ( 1.5 - (x_2 * result * result )) # Newton's method (again)
        return result
    end
    function fast_inv_sqrt(x::Float32)
        x_2 = 0.5f0 * x
        result = Base.sub_int(1.321202f19, Base.lshr_int(x,1)) # Floating point magic
        result *= ( 1.5f0 - (x_2 * result * result) ) # Newton's method
        result *= ( 1.5f0 - (x_2 * result * result) ) # Newton's method (again)
        return result
    end
    export fast_inv_sqrt

## --- Remove non-positive numbers

    function positive!(a::DenseArray{<:AbstractFloat})
        @inbounds for i in eachindex(a)
            if !(a[i] > 0)
                a[i] = NaN
            end
        end
        return a
    end
    export positive!

## --- Rescale an AbstractArray between a new minimum and maximum

    """
    ```julia
    rescale(y, min::Number=0, max::Number=1)
    ```
    Rescale a collection of numbers `y` between a new minimum `min` and new maximum `max`

    ### Examples
    ```julia
    julia> rescale(1:5)
    5-element Vector{Float64}:
    0.0
    0.25
    0.5
    0.75
    1.0

    julia> rescale(1:5, -1, 0)
    5-element Vector{Float64}:
    -1.0
    -0.75
    -0.5
    -0.25
    0.0
    ```
    """
    function rescale(y::AbstractArray, min::Number=0, max::Number=1)
        obsmin = nanminimum(y)
        y = float.(y) .- obsmin
        obsmax = nanmaximum(y)
        y ./= obsmax
        y .+= min
        y .*= (max-min)
    end
    function rescale(y::AbstractRange, min::Number=0, max::Number=1)
        obsmin = minimum(y)
        y = float(y) .- obsmin
        obsmax = maximum(y)
        y = y ./ obsmax
        y = y .+ min
        y = y .* (max-min)
    end
    function rescale(y, min::Number=0, max::Number=1)
        obsmin = nanminimum(y)
        y = float.(y) .- obsmin
        obsmax = nanmaximum(y)
        y = y ./ obsmax
        y = y .+ min
        y = y .* (max-min)
    end
    export rescale

## --- Some mathematical constants

    const SQRT2 = sqrt(2)
    const SQRT2PI = sqrt(2*pi)
    const INVSQRT2 = 1/sqrt(2)
    const AN = Union{AbstractArray{<:Number},Number}

## --- Gaussian distribution functions

    """
    ```julia
    normpdf(mu,sigma,x)
    ```
    Probability density function of the Normal (Gaussian) distribution

    ``ℯ^{-(x-μ)^2 / (2σ^2)} / σ√2π``

    with mean `mu` and standard deviation `sigma`, evaluated at `x`
    """
    @inline normpdf(mu,sigma,x) = exp(-(x-mu)*(x-mu) / (2*sigma*sigma)) / (SQRT2PI*sigma)
    @inline normpdf(mu::Number,sigma::Number,x::Number) = exp(-(x-mu)*(x-mu) / (2*sigma*sigma)) / (SQRT2PI*sigma)
    normpdf(mu::AN,sigma::AN,x::AN) = @fastmath @. exp(-(x-mu)*(x-mu) / (2*sigma*sigma)) / (SQRT2PI*sigma)
    export normpdf

    """
    ```julia
    normpdf_ll(mu, sigma, x)
    ```
    Fast log likelihood proportional to the natural logarithm of the probability
    density function of a Normal (Gaussian) distribution with mean `mu` and
    standard deviation `sigma`, evaluated at `x`.

    If `x`, [`mu`, and `sigma`] are given as arrays, the sum of the log likelihood
    over all `x` will be returned.

    See also `normpdf`, `normlogpdf`
    """
    @inline normpdf_ll(mu,sigma,x) = -(x-mu)*(x-mu) / (2*sigma*sigma)
    function normpdf_ll(mu::Number,sigma::Number,x::AbstractArray)
        inv_s2 = 1/(2*sigma*sigma)
        ll = zero(typeof(inv_s2))
        @inbounds @fastmath @simd ivdep for i ∈ eachindex(x)
            ll -= (x[i]-mu)*(x[i]-mu) * inv_s2
        end
        return ll
    end
    function normpdf_ll(mu::AbstractArray,sigma::Number,x::AbstractArray)
        inv_s2 = 1/(2*sigma*sigma)
        ll = zero(typeof(inv_s2))
        @inbounds @fastmath @simd ivdep for i ∈ eachindex(x, mu)
            ll -= (x[i]-mu[i])*(x[i]-mu[i]) * inv_s2
        end
        return ll
    end
    function normpdf_ll(mu::Number,sigma::AbstractArray,x::AbstractArray)
        ll = zero(float(eltype(sigma)))
        @inbounds @fastmath @simd ivdep for i ∈ eachindex(x, sigma)
            ll -= (x[i]-mu)*(x[i]-mu) / (2*sigma[i]*sigma[i])
        end
        return ll
    end
    function normpdf_ll(mu::AbstractArray,sigma::AbstractArray,x::AbstractArray)
        ll = zero(float(eltype(sigma)))
        @inbounds @fastmath @simd ivdep for i ∈ eachindex(x, mu, sigma)
            ll -= (x[i]-mu[i])*(x[i]-mu[i]) / (2*sigma[i]*sigma[i])
        end
        return ll
    end
    export normpdf_ll


    """
    ```julia
    normlogpdf(mu, sigma, x)
    ```
    The natural logarithm of the probability density function of a Normal (Gaussian) 
    distribution with mean `mu` and standard deviation `sigma`, evaluated at `x`.

    If `x`, [`mu`, and `sigma`] are given as arrays, the sum of the log probability density
    over all `x` will be returned.

    See also `normpdf`, `normlogpdf`
    """
    @inline normlogpdf(mu,sigma,x) = -(x-mu)*(x-mu) / (2*sigma*sigma) - log(SQRT2PI*sigma)
    function normlogpdf(mu::Number,sigma::Number,x::AbstractArray)
        inv_s2 = 1/(2*sigma*sigma)
        ll = zero(typeof(inv_s2))
        @inbounds @fastmath @simd ivdep for i ∈ eachindex(x)
            ll -= (x[i]-mu)*(x[i]-mu) * inv_s2
        end
        return ll - length(x)*log(SQRT2PI*sigma)
    end
    function normlogpdf(mu::AbstractArray,sigma::Number,x::AbstractArray)
        inv_s2 = 1/(2*sigma*sigma)
        ll = zero(typeof(inv_s2))
        @inbounds @fastmath @simd ivdep for i ∈ eachindex(x, mu)
            ll -= (x[i]-mu[i])*(x[i]-mu[i]) * inv_s2
        end
        return ll - log(SQRT2PI*sigma)*length(x)
    end
    function normlogpdf(mu::Number,sigma::AbstractArray,x::AbstractArray)
        ll = zero(float(eltype(sigma)))
        @inbounds @fastmath @simd ivdep for i ∈ eachindex(x, sigma)
            ll -= (x[i]-mu)*(x[i]-mu) / (2*sigma[i]*sigma[i]) + log(SQRT2PI*sigma[i])
        end
        return ll
    end
    function normlogpdf(mu::AbstractArray,sigma::AbstractArray,x::AbstractArray)
        ll = zero(float(eltype(sigma)))
        @inbounds @fastmath @simd ivdep for i ∈ eachindex(x, mu, sigma)
            ll -= (x[i]-mu[i])*(x[i]-mu[i]) / (2*sigma[i]*sigma[i]) + log(SQRT2PI*sigma[i])
        end
        return ll
    end
    export normlogpdf

    """
    ```julia
    normcdf(mu,sigma,x)
    ```
    Cumulative distribution function of the Normal (Gaussian) distribution

    ``1/2 + erf(\frac{x-μ}{σ√2})/2``

    with mean `mu` and standard deviation `sigma`, evaluated at `x`.
    """
    @inline normcdf(mu,sigma,x) = 0.5 + 0.5 * erf((x-mu) / (sigma*SQRT2))
    @inline normcdf(mu::Number,sigma::Number,x::Number) = 0.5 + 0.5 * erf((x-mu) / (sigma*SQRT2))
    normcdf(mu::AN,sigma::AN,x::AN) = @fastmath @. 0.5 + 0.5 * erf((x-mu) / (sigma*SQRT2))
    export normcdf


    """
    ```julia
    normcdf_ll(mu, sigma, x)
    ```
    Fast log likelihood proportional to the natural logarithm of the cumulative
    distribution function of a Normal (Gaussian) distribution with mean `mu` and
    standard deviation `sigma`, evaluated at `x`.

    If `x`, [`mu`, and `sigma`] are given as arrays, the sum of the log likelihood
    over all `x` will be returned.

    See also `normcdf`
    """
    @inline function normcdf_ll(xₛ::Number)
        if xₛ < -1.0
            return log(0.5*erfcx(-xₛ * INVSQRT2)) - 0.5*abs2(xₛ)
        else
            return log1p(-0.5*erfc(xₛ * INVSQRT2))
        end
    end
    function normcdf_ll(xₛ::AbstractArray)
        ll = zero(float(eltype(xₛ)))
        @inbounds for i ∈ eachindex(xₛ)
            ll += normcdf_ll(xₛ[i])
        end
        return ll
    end
    @inline function normcdf_ll(mu::Number, sigma::Number, x::Number)
        xₛ = (x - mu) / sigma
        return normcdf_ll(xₛ)
    end
    function normcdf_ll(mu::Number,sigma::Number,x::AbstractArray)
        inv_sigma = 1/sigma
        ll = zero(typeof(inv_sigma))
        @inbounds for i ∈ eachindex(x)
            xₛ = (x[i] - mu) * inv_sigma
            ll += normcdf_ll(xₛ)
        end
        return ll
    end
    function normcdf_ll(mu::AbstractArray,sigma::Number,x::AbstractArray)
        inv_sigma = 1/sigma
        ll = zero(typeof(inv_sigma))
        @inbounds for i ∈ eachindex(x)
            xₛ = (x[i] - mu[i]) * inv_sigma
            ll += normcdf_ll(xₛ)
        end
        return ll
    end
    function normcdf_ll(mu::Number,sigma::AbstractArray,x::AbstractArray)
        ll = zero(float(eltype(sigma)))
        @inbounds for i ∈ eachindex(x)
            xₛ = (x[i] - mu) / sigma[i]
            ll += normcdf_ll(xₛ)
        end
        return ll
    end
    function normcdf_ll(mu::AbstractArray,sigma::AbstractArray,x::AbstractArray)
        ll = zero(float(eltype(sigma)))
        @inbounds for i ∈ eachindex(x)
            xₛ = (x[i] - mu[i]) / sigma[i]
            ll += normcdf_ll(xₛ)
        end
        return ll
    end
    export normcdf_ll


    """
    ```julia
    normcdf_ll!(mu, sigma, x)
    ```
    Fast log likelihood proportional to the natural logarithm of the cumulative
    distribution function of a Normal (Gaussian) distribution with mean `mu` and
    standard deviation `sigma`, evaluated at `x`.

    As `normcdf_ll`, but in-place (using `x` as a buffer).
    """
    function normcdf_ll!(xₛ::AbstractArray)
        @inbounds for i ∈ eachindex(xₛ)
            xₛ[i] = normcdf_ll(xₛ[i])
        end
        ll = zero(float(eltype(xₛ)))
        @inbounds @fastmath for i ∈ eachindex(xₛ)
            ll += xₛ[i]
        end
        return ll
    end
    function normcdf_ll!(mu::Number,sigma::Number,x::AbstractArray)
        inv_sigma = 1/sigma
        @inbounds for i ∈ eachindex(x)
            xₛ = (x[i] - mu) * inv_sigma
            x[i] = normcdf_ll(xₛ)
        end
        ll = zero(typeof(inv_sigma))
        @inbounds @fastmath for i ∈ eachindex(x)
            ll += x[i]
        end
        return ll
    end
    function normcdf_ll!(mu::AbstractArray,sigma::Number,x::AbstractArray)
        inv_sigma = 1/sigma
        @inbounds for i ∈ eachindex(x)
            xₛ = (x[i] - mu[i]) * inv_sigma
            x[i] = normcdf_ll(xₛ)
        end
        ll = zero(typeof(inv_sigma))
        @inbounds @fastmath for i ∈ eachindex(x)
            ll += x[i]
        end
        return ll
    end
    function normcdf_ll!(mu::Number,sigma::AbstractArray,x::AbstractArray)
        @inbounds for i ∈ eachindex(x)
            xₛ = (x[i] - mu) / sigma[i]
            x[i] = normcdf_ll(xₛ)
        end
        ll = zero(float(eltype(sigma)))
        @inbounds @fastmath for i ∈ eachindex(x)
            ll += x[i]
        end
        return ll
    end
    function normcdf_ll!(mu::AbstractArray,sigma::AbstractArray,x::AbstractArray)
        @inbounds for i ∈ eachindex(x)
            xₛ = (x[i] - mu[i]) / sigma[i]
            x[i] = normcdf_ll(xₛ)
        end
        ll = zero(float(eltype(sigma)))
        @inbounds @fastmath for i ∈ eachindex(x)
            ll += x[i]
        end
        return ll
    end
    export normcdf_ll!

    """
    ```julia
    normcdf!(result,mu,sigma,x)
    ```
    In-place version of `normcdf`
    """
    function normcdf!(result::DenseArray, mu::Number, sigma::Number, x::AbstractArray)
        T = eltype(result)
        inv_sigma_sqrt2 = one(T)/(sigma*T(SQRT2))
        @inbounds @fastmath for i ∈ eachindex(x,result)
            result[i] = T(0.5) + T(0.5) * erf((x[i]-mu) * inv_sigma_sqrt2)
        end
        return result
    end
    export normcdf!

    """
    ```julia
    norm_quantile(F::Number)
    ```
    How far away from the mean (in units of sigma) should we expect proportion
    F of the samples to fall in a standard Gaussian (Normal[0,1]) distribution
    """
    @inline norm_quantile(F) = SQRT2*erfinv(2*F-1)
    export norm_quantile

    """
    ```julia
    norm_width(N::Number)
    ```
    How dispersed (in units of sigma) should we expect a sample of N numbers
    drawn from a standard Gaussian (Normal[0,1]) distribution to be?
    """
    @inline norm_width(N) = 2*norm_quantile(1 - 1/(2N))
    export norm_width

    """
    ```julia
    normproduct(μ1, σ1, μ2, σ2)
    ```
    The integral of the product of two normal distributions N[μ1,σ1] * N[μ2,σ2].
    This is itself just another Normal distribution! Specifically, one with
    variance σ1^2 + σ2^2, evaluated at distance |μ1-μ2| from the mean
    """
    normproduct(μ1, σ1, μ2, σ2) = normpdf(μ1, sqrt.(σ1.*σ1 + σ2.*σ2), μ2)
    export normproduct

    """
    ```julia
    normproduct_ll(μ1, σ1, μ2, σ2)
    ```
    Fast log likelihood proportional to the integral of N[μ1,σ1] * N[μ2,σ2]
    As `normlogproduct`, but using the fast log likelihood of a Normal distribution
    (i.e., without the preexponential terms).
    """
    normproduct_ll(μ1, σ1, μ2, σ2) = normpdf_ll(μ1, sqrt.(σ1.*σ1 + σ2.*σ2), μ2)
    export normproduct_ll

    """
    ```julia
    normlogproduct(μ1, σ1, μ2, σ2)
    ```
    The logarithm of the integral of N[μ1,σ1] * N[μ2,σ2]
    """
    normlogproduct(μ1, σ1, μ2, σ2) = normlogpdf(μ1, sqrt.(σ1.*σ1 + σ2.*σ2), μ2)
    export normlogproduct


## --- Geometry

    """
    ```julia
    inpolygon(x,y,point)
    ```
    Check if a 2D polygon defined by the arrays `x`, `y` contains a given `point`.
    Returns boolean (true or false)

    ### Examples
    ```julia
    julia> x = [0, 1, 1, 0];

    julia> y = [0, 0, 1, 1];

    julia> inpolygon(x, y, (0.5,0.5))
    true

    julia> inpolygon(x, y, (0.5,1.5))
    false
    ```
    """
    function inpolygon(x,y,point)
        # Check that we have the right kind of input data
        if length(x) != length(y)
            error("polygon must have equal number of x and y points\n")
        end
        if length(x) < 3
            error("polygon must have at least 3 points\n")
        end
        if length(point) != 2
            error("point must be an ordered pair (x,y)\n")
        end

        # Extract x and y data of point
        point_x = point[1]
        point_y = point[2]
        # For first point, previous point is last
        x_here = x[end]
        y_here = y[end]
        # Ensure we are not sitting parallel to a vertex by infinitessimally moving the point
        if y_here == point_y
            point_y = nextfloat(float(point_y))
        end
        if x_here == point_x
            point_x = nextfloat(float(point_x))
        end

        # Check how many times a line projected right along x-axis from point intersects the polygon
        intersections = 0
        @inbounds for i ∈ eachindex(x)
            # Recycle our vertex
            x_last = copy(x_here)
            y_last = copy(y_here)

            # Get a new vertex
            x_here = x[i]
            y_here = y[i]

            # Ensure we are not sitting parallel to a vertex by infinitessimally moving the point
            if y_here == point_y
                point_y = nextfloat(float(point_y))
            end
            if x_here == point_x
                point_x = nextfloat(float(point_x))
            end

            if y_last > point_y && y_here > point_y
                # If both ys above point, no intersection
                continue
            elseif y_last < point_y && y_here < point_y
                # If both ys below point, no intersection
                continue
            elseif x_last < point_x && x_here < point_x
                # If both x's left of point, no intersection
                continue
            elseif x_last > point_x && x_here > point_x
                # By elimination
                # We have one y above and y below our point
                # If both y's are right of line, then definite intersection
                intersections += 1
                continue
            else
                # By elimination
                # One y above and one y below
                # One x to the right and one x to the left
                # We must project
                dy = y_here - y_last
                if abs(dy) > 0
                    dx = x_here - x_last
                    inv_slope = dx / dy
                    x_proj = x_last + (point_y - y_last) * inv_slope
                    if x_proj > point_x
                        intersections += 1
                    end
                end
            end
        end

        # If number of intersections is odd, point is in the polygon
        return Bool(mod(intersections,2))
    end
    export inpolygon


    """
    ```julia
    (columns, rows) = find_grid_inpolygon(grid_x, grid_y, poly_x, poly_y)
    ```
    Find the indexes of grid points that fall within a polygon for a grid with
    cell centers given by grid_x (j-columns of grid) and grid_y (i-rows of grid).
    Returns a list of rows and columns in the polygon

    ### Examples
    ```julia
    julia> grid_x = -1.5:1/3:1.5;

    julia> grid_y = -1.5:1/3:1.5;

    julia> cols,rows = find_grid_inpolygon(gridx, gridy, [-.4,.4,.4,-.4],[.4,.4,-.4,-.4])
    ([5, 5, 6, 6], [5, 6, 5, 6])

    julia> grid_x[cols]
    4-element Vector{Float64}:
     -0.16666666666666666
     -0.16666666666666666
      0.16666666666666666
      0.16666666666666666

    julia> grid_y[rows]
    4-element Vector{Float64}:
     -0.16666666666666666
      0.16666666666666666
     -0.16666666666666666
      0.16666666666666666
    """
    function find_grid_inpolygon(grid_x, grid_y, poly_x, poly_y)
        # Check that we have the right kind of input data
        if length(poly_x) != length(poly_y)
            error("polygon must have equal number of x and y points\n")
        end
        if length(poly_x) < 3
            error("polygon must have at least 3 points\n")
        end

        # Find maximum x and y range of polygon
        (xmin, xmax) = extrema(poly_x)
        (ymin, ymax) = extrema(poly_y)

        # Find the matrix indices within the range of the polygon (if any)
        column_inrange = findall((grid_x .>= xmin) .& (grid_x .<= xmax))
        row_inrange = findall((grid_y .>= ymin) .& (grid_y .<= ymax))

        # Keep a list of matrix indexes in the polygon
        row = Array{Int}(undef,length(column_inrange) * length(row_inrange))
        column = Array{Int}(undef,length(column_inrange) * length(row_inrange))
        n = 0
        for j ∈ eachindex(column_inrange)
            for i ∈ eachindex(row_inrange)
                point = (grid_x[column_inrange[j]], grid_y[row_inrange[i]])
                if inpolygon(poly_x, poly_y, point)
                    n += 1
                    row[n] = row_inrange[i]
                    column[n] = column_inrange[j]
                end
            end
        end

        return (column[1:n], row[1:n])
    end
    export find_grid_inpolygon


    """
    ```julia
    arcdistance(latᵢ,lonᵢ,lat,lon)
    ```
    Calculate the distance on a sphere between the point (`latᵢ`,`lonᵢ`) and any
    number of points in (`lat`,`lon`).
    Latitude and Longitude should be specified in decimal degrees
    """
    function arcdistance(latᵢ,lonᵢ,lat,lon)
        @assert eachindex(latᵢ) == eachindex(lonᵢ)
        @assert eachindex(lat) == eachindex(lon)

        # Argument for acos()
        arg = @. sin(latᵢ * pi/180) * sin(lat * pi/180) + cos(latᵢ*pi/180) * cos(lat * pi/180)*cos((lonᵢ - lon) * pi/180)

        # Avoid domain errors from imprecise sine and cosine math
        @inbounds for i in eachindex(arg)
            if arg[i] < -1
                arg[i] = -1
            elseif arg[i] > 1
                arg[i] = 1
            end
        end

        # Calculate angular distance
        theta = 180/pi .* acos.(arg)
        return theta
    end
    export arcdistance

    """
    ```julia
    minarcdistance(latᵢ,lonᵢ,lat,lon)
    ```
    Return the smallest non-`NaN` arcdistance (i.e. distance on a sphere in arc degrees)
    between a given point (`latᵢ[i]`,`lonᵢ[i]`) and any point in (`lat`,`lon`)
    for each `i` in `eachindex(latᵢ, lonᵢ)`.

    Latitude and Longitude should be specified in decimal degrees
    """
    function minarcdistance(latᵢ,lonᵢ,lat,lon)
        @assert eachindex(latᵢ) == eachindex(lonᵢ)
        @assert eachindex(lat) == eachindex(lon)

        # Precalculate some shared factors
        sli = sin.(latᵢ .* pi/180)
        sl = sin.(lat .* pi/180)
        cli = cos.(latᵢ*pi/180)
        cl = cos.(lat .* pi/180)

        thetamin = fill(NaN, size(latᵢ))
        @inbounds for i in eachindex(latᵢ)
            for j in eachindex(lon)
                arg = sli[i] * sl[j] + cli[i] * cl[j] * cos((lonᵢ[i] - lon[j]) * pi/180)
                if arg < -1
                    arg = -1.0
                elseif arg > 1
                    arg = 1.0
                end
                θᵢⱼ = 180/pi * acos(arg)
                if !(θᵢⱼ >= thetamin[i])
                    thetamin[i] = θᵢⱼ
                end
            end
        end
        return thetamin
    end
    export minarcdistance


## --- Linear regression

    """
    ```julia
    (a,b) = linreg(x::AbstractVector, y::AbstractVector)
    ```
    Returns the coefficients for a simple linear least-squares regression of
    the form `y = a + bx`

    ### Examples
    ```
    julia> a, b = linreg(1:10, 1:10)
    2-element Vector{Float64}:
     -1.19542133983862e-15
      1.0

    julia> isapprox(a, 0, atol = 1e-12)
    true

    julia> isapprox(b, 1, atol = 1e-12)
    true
    ```
    """
    function linreg(x::AbstractVector{T}, y::AbstractVector{<:Number}) where {T<:Number}
        A = similar(x, length(x), 2)
        A[:,1] .= one(T)
        A[:,2] .= x
        return A\y
    end
    export linreg


## --- End of File

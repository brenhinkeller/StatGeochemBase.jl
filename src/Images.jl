## --- Calculate a 2d path density image from x-y path distributions

    """
    ```julia
    image_from_paths(xpoints::AbstractMatrix, ypoints::AbstractMatrix; 
        xresolution::Int=1800, 
        yresolution::Int=1200, 
        xrange=nanextrema(xpoints), 
        yrange=nanextrema(ypoints),
    )
    ```
    Produce a 2d image (histogram) of path densities given a set of x-y points 
    stored column-wise in separate matrices `xpoints` and `ypoints`.
    `x` is assumed to be the independent variable (i.e., the x-y points are
    plotted in order of increasing `x`).

    See also `image_from_paths!`

    ### Examples
    ```julia
    julia> nsims = 1000
    1000

    julia> xdist = rand(10, nsims);

    julia> ydist = rand(10, nsims);

    julia> imgcounts, xbincenters, ybincenters = image_from_paths(xpaths, ypaths; xresolution=50, yresolution=50)
    ([8 15 … 9 11; 12 4 … 14 11; … ; 9 15 … 9 9; 10 17 … 10 14], -3.715247394908744:0.1523101612461508:3.747950506152645, -3.86556932981587:0.1497772964483582:3.4735181961536803)
    ```
    """
    image_from_paths(xpoints, ypoints; kwargs...) = image_from_paths!(copy(xpoints), copy(ypoints); kwargs...)
    export image_from_paths

    """
    ```julia
    image_from_paths!([imgcounts::AbstractMatrix], xpoints::AbstractMatrix, ypoints::AbstractMatrix; 
        xresolution::Int=1800, 
        yresolution::Int=1200, 
        xrange=nanextrema(xpoints), 
        yrange=nanextrema(ypoints),
        method=:nearest,
    )
    ```
    As `image_from_paths`, but may sort `xpoints`  (as required for interpolation) 
    in-place rather than making a copy.

    Optionally, a result matrix `imgcounts` may be supplied for fully in-place operation, 
    in which case `yresolution, xresolution = size(imgcounts)`

    Available line-drawing `method`s include 
        `:interpolate`, which interpolates the independent variable (y) as a function of the dependent variable (x)
        `:closest`, which attempts to find the 
    """
    image_from_paths!(xpoints::AbstractMatrix, ypoints::AbstractMatrix; xresolution::Int=1800, yresolution::Int=1200, kwargs...) = image_from_paths!(zeros(Int, yresolution, xresolution), xpoints, ypoints; kwargs...)
    function image_from_paths!(imgcounts::AbstractMatrix, xpoints::AbstractMatrix, ypoints::AbstractMatrix; 
            xrange=nanextrema(xpoints), 
            yrange=nanextrema(ypoints),
            method::Symbol=:nearest,
        )

        # Bin edges and centers
        yresolution, xresolution = size(imgcounts)
        xbinedges = range(first(xrange), last(xrange), length=xresolution+1)
        xq = cntr(xbinedges)
        ybinedges = range(first(yrange), last(yrange), length=yresolution+1)
        yq = cntr(ybinedges)
        # Decide method
        (method === :interpolate || method === :nearest) || @warn "Method $method not recognized, falling back to :nearest"
        if method===:interpolate
            interpolate_pixels!(imgcounts, xpoints, ypoints, xq, ybinedges)
        else
            nearest_pixels!(imgcounts, xpoints, ypoints, xq, yq)
        end
        return imgcounts, xq, yq
    end
    export image_from_paths!

    function interpolate_pixels!(imgcounts, xpoints, ypoints, xq, ybinedges)
        @assert axes(xpoints, 1) == axes(ypoints, 1)
        nsims = size(xpoints, 2)
        @assert axes(xpoints, 2) == axes(ypoints, 2) == Base.OneTo(nsims)
        yresolution, xresolution = size(imgcounts)
        @assert length(xq) == xresolution
        @assert length(ybinedges) == yresolution + 1

        # Use batches to avoid allocating massive intermediate arrays
        batchsize = min(nsims, 1000)
        interpydist = zeros(xresolution, batchsize)
        interpyknots = zeros(Int, xresolution)

        # Loop through batches
        n₀ = 0
        while n₀+batchsize <= nsims
            # Interpolate paths to match image resolution
            for nᵢ in Base.OneTo(batchsize)
                n = nᵢ+n₀
                linterp1s!(view(interpydist, :, nᵢ), interpyknots, view(xpoints, :, n), view(ypoints, :, n), xq)
            end

            # Calculate composite image
            for i in Base.OneTo(xresolution) # scan by x (one column at a time)
                histcounts!(view(imgcounts,:,i), view(interpydist,i,:), ybinedges)
            end
            n₀ += batchsize
        end
        nextra = nsims-n₀
        if nextra>0
            # Interpolate paths to match image resolution
            for nᵢ in Base.OneTo(nextra)
                n = nᵢ+n₀
                linterp1s!(view(interpydist, :, nᵢ), interpyknots, view(xpoints, :, n), view(ypoints, :, n), xq)
            end

            # Calculate composite image
            for i in Base.OneTo(xresolution) # scan by x (one column at a time)
                histcounts!(view(imgcounts,:,i), view(interpydist,i,1:nextra), ybinedges)
            end
        end
        return imgcounts
    end
    function nearest_pixels!(imgcounts, xpoints, ypoints, xq, yq)
        npoints = size(xpoints, 1)
        @assert axes(xpoints, 1) == axes(ypoints, 1) == Base.OneTo(npoints)
        nsims = size(xpoints, 2)
        @assert axes(xpoints, 2) == axes(ypoints, 2) == Base.OneTo(nsims)
        @assert eachindex(yq) == axes(imgcounts,1)
        @assert eachindex(xq) == axes(imgcounts,2)
        
        # Allocate buffers
        x = zeros(size(xpoints,1))
        y = zeros(size(xpoints,1))

        # Loop through each set of paths
        for n in Base.OneTo(nsims)
            # Copy x and y points to buffers
            copyto!(x, 1, xpoints, 1+(n-1)*npoints, npoints)
            copyto!(y, 1, ypoints, 1+(n-1)*npoints, npoints)
            # Sort x and permute y to match, in-place
            nanargsort!(y, x) 
            # Convert to image coordinates
            x .-= first(xq)
            x ./= step(xq)
            x .+= firstindex(xq)
            y .-= first(yq)
            y ./= step(yq)
            y .+= firstindex(yq)
            for i in Base.OneTo(npoints-1)
                Δx = x[i+1] - x[i]
                isnan(Δx) && break
                Δy = y[i+1] - y[i]
                isnan(Δy) && break
                s = Δy/Δx
                if abs(s) > 1
                    # Steep line
                    # Raster vertically (by y) and find closest x pixel at each y
                    sinv = Δx/Δy
                    p_start, p_end = (s > 1) ? (i, i+1) : (i+1, i)
                    y_start, y_end = y[p_start], y[p_end]
                    xi = x[p_start]
                    iy = ceil(Int, y_start)
                    while iy < y_end
                        ix = round(Int, xi)
                        if (firstindex(xq) <= ix <= lastindex(xq)) && (firstindex(yq) <= iy <= lastindex(yq))
                            imgcounts[iy, ix] += 1
                        end
                        iy += 1
                        xi += sinv
                    end
                else
                    # Shallow line
                    # Raster horizontally (by x) and find closest y pixel at each x
                    x_start, x_end = x[i], x[i+1]
                    yi = y[i]
                    ix = ceil(Int, x_start)
                    while ix < x_end
                        iy = round(Int, yi)
                        if (firstindex(xq) <= ix <= lastindex(xq)) && (firstindex(yq) <= iy <= lastindex(yq))
                            imgcounts[iy, ix] += 1
                        end
                        ix += 1
                        yi += s
                    end
                end
            end
        end

        return imgcounts
    end

## --- Map colormaps to images

    """
    ```julia
    imsc(A::AbstractArray, colormap::AbstractVector=viridis, cmin=nanminimum(A), cmax=nanmaximum(A))
    ```
    Convert a matrix `A` to an image (an array of Colors.jl colors) using the
    specified `colormap` (default `viridis`), optionally scaled between `cmin`
    and `cmax`.

    ### Examples
    ```julia
    julia> A = rand(3,3)
    3×3 Matrix{Float64}:
     0.39147   0.553489  0.351628
     0.331786  0.343836  0.824674
     0.639233  0.558113  0.965627

    julia> img = imsc(A) # N.B. will display as image if `using ImageInTerminal`
    3×3 Array{RGB{N0f8},2} with eltype ColorTypes.RGB{FixedPointNumbers.N0f8}:
     RGB{N0f8}(0.282,0.137,0.455)  …  RGB{N0f8}(0.278,0.051,0.376)
     RGB{N0f8}(0.267,0.004,0.329)     RGB{N0f8}(0.431,0.808,0.345)
     RGB{N0f8}(0.133,0.553,0.553)     RGB{N0f8}(0.992,0.906,0.145)

    julia> using Images; save("img.png", img) # Save to file as PNG

    julia> using Plots; plot(0:3, 0:3, img) # Plot with specified x and y axes
    ```
    """
    function imsc(A::AbstractArray, colormap::AbstractVector=viridis, cmin=nanminimum(A), cmax=nanmaximum(A))
        Nc = length(colormap)
        crange = cmax - cmin
        return  A .|> x -> colormap[isnan(x) ? 1 : ceil(UInt, min(max(Nc*(x-cmin)/crange, 1), Nc))]
    end
    export imsc

    """
    ```julia
    imsci(A::AbstractArray, colormap::AbstractVector=viridis, cmin=nanminimum(A), cmax=nanmaximum(A))
    ```
    Convert a matrix `A` to an indirect array image (an IndirectArray of Colors.jl
    colors) using the specified `colormap` (default `viridis`), optionally scaled
    between `cmin` and `cmax`.

    As `imsc`, but returns an IndirectArray; slightly more space efficient for
    small colormaps, but with computational cost.

    ### Examples
    ```julia
    julia> A = rand(3,3)
    3×3 Matrix{Float64}:
     0.39147   0.553489  0.351628
     0.331786  0.343836  0.824674
     0.639233  0.558113  0.965627

     julia> img = imsci(A)
     3×3 IndirectArrays.IndirectArray{RGB{N0f8}, 2, UInt64, Matrix{UInt64}, Vector{RGB{N0f8}}}:
      RGB{N0f8}(0.282,0.137,0.455)  …  RGB{N0f8}(0.278,0.051,0.376)
      RGB{N0f8}(0.267,0.004,0.329)     RGB{N0f8}(0.431,0.808,0.345)
      RGB{N0f8}(0.133,0.553,0.553)     RGB{N0f8}(0.992,0.906,0.145)

    julia> using Images; save("img.png", img) # Save to file as PNG

    julia> using Plots; plot(0:3, 0:3, img) # Plot with specified x and y axes
    ```
    """
    function imsci(A::AbstractArray,colormap::AbstractArray=viridis,cmin::Number=nanminimum(A),cmax::Number=nanmaximum(A))
        Nc = length(colormap)
        crange = cmax - cmin
        return IndirectArray(A .|> x -> isnan(x) ? 1 : ceil(UInt, min(max(Nc*(x-cmin)/crange, 1), Nc)), colormap)
    end
    export imsci


## -- End of File

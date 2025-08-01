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
    )
    ```
    As `image_from_paths`, but will sort `xpoints`  (as required for interpolation) 
    in-place rather than making a copy.

    Optionally, a result matrix `imgcounts` may be supplied for fully in-place operation, 
    in which case `yresolution, xresolution = size(imgcounts)`
    """
    image_from_paths!(xpoints::AbstractMatrix, ypoints::AbstractMatrix; xresolution::Int=1800, yresolution::Int=1200, kwargs...) = image_from_paths!(zeros(Int, yresolution, xresolution), xpoints, ypoints; kwargs...)
    function image_from_paths!(imgcounts::AbstractMatrix, xpoints::AbstractMatrix, ypoints::AbstractMatrix; xrange=nanextrema(xpoints), yrange=nanextrema(ypoints))
        yresolution, xresolution = size(imgcounts)
        @assert axes(xpoints, 1) == axes(ypoints, 1)
        nsims = size(xpoints, 2)
        @assert axes(xpoints, 2) == axes(ypoints, 2) == Base.OneTo(nsims)

        # Bin edges and centers
        xbinedges = range(first(xrange), last(xrange), length=xresolution+1)
        xq = cntr(xbinedges)
        ybinedges = range(first(yrange), last(yrange), length=yresolution+1)
        yq = cntr(ybinedges)

        # Use batches to avoid allocating massive intermediate arrays
        batchsize = min(nsims, 1000)
        interpydist = zeros(xresolution, batchsize)
        interpyknots = zeros(Int, xresolution)
        # interpxdist = zeros(yresolution, batchsize)
        # interpxknots = zeros(Int, yresolution)

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

        return imgcounts, xq, yq
    end
    export image_from_paths!

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

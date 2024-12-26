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
    Produce a 2d image (histogram) of path densities given a set of x and y paths 
    stored column-wise in separate matrices `xpoints` and `ypoints`.

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
    image_from_paths(xpoints, ypoints; xresolution::Int=1800, yresolution::Int=1200, xrange=nanextrema(xpoints), yrange=nanextrema(ypoints)) = image_from_paths!(copy(xpoints), copy(ypoints); xresolution, yresolution, xrange, yrange)
    export image_from_paths

    function image_from_paths!(xpoints::AbstractMatrix, ypoints::AbstractMatrix; xresolution::Int=1800, yresolution::Int=1200, xrange=nanextrema(xpoints), yrange=nanextrema(ypoints))
        @assert axes(xpoints, 1) == axes(ypoints, 1)
        nsims = size(xpoints, 2)
        @assert axes(xpoints, 2) == axes(ypoints, 2) == Base.OneTo(nsims)

        # Interpolate paths to match image resolution
        xbinedges = range(first(xrange), last(xrange), length=xresolution+1)
        xq = cntr(xbinedges)
        ybinedges = range(first(yrange), last(yrange), length=yresolution+1)
        yq = cntr(ybinedges)
        xinterpdist = Array{Float64}(undef, xresolution, nsims)
        for i in Base.OneTo(nsims)
            linterp1s!(view(xinterpdist,:,i), view(xpoints, :, i), view(ypoints, :, i), xq)
        end
        yinterpdist = Array{Float64}(undef, yresolution, nsims)
        for i in Base.OneTo(nsims)
            linterp1s!(view(yinterpdist,:,i), view(xpoints, :, i), view(ypoints, :, i), yq)
        end

        # Calculate composite image, scanning both by x and y
        imgcounts = zeros(Int, yresolution, xresolution)
        for i in Base.OneTo(xresolution) # scan by x (one column at a time)
            histcounts!(view(imgcounts,:,i), view(xinterpdist,i,:), ybinedges)
        end
        for j in Base.OneTo(yresolution) # scan by y (one row at a time)
            histcounts!(view(imgcounts,j,:), view(yinterpdist,j,:), xbinedges)
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

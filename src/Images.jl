## --- Map colormaps to images

    # Convert matrix to image using colormap
    function imsc(matrix::AbstractArray,colormap::AbstractArray=viridis,cmin::Number=0,cmax::Number=0)
        Nc = length(colormap)
        if cmin>=cmax
            cmin = nanminimum(matrix)
            cmax = nanmaximum(matrix)
        end
        crange = cmax - cmin
        return  matrix .|> x -> colormap[isnan(x) ? 1 : ceil(UInt, min(max(Nc*(x-cmin)/crange, 1), Nc))]
    end
    export imsc

    # Convert matrix to indirect array image using colormap
    function imsci(matrix::AbstractArray,colormap::AbstractArray=viridis,cmin::Number=0,cmax::Number=0)
        Nc = length(colormap)
        if cmin>=cmax
            cmin = nanminimum(matrix)
            cmax = nanmaximum(matrix)
        end
        crange = cmax - cmin
        return IndirectArray(matrix .|> x -> isnan(x) ? 1 : ceil(UInt, min(max(Nc*(x-cmin)/crange, 1), Nc)), colormap)
    end
    export imsci


## -- End of File

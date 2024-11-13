module StatGeochemBase

    using NaNStatistics

    using SpecialFunctions: erf, erfc, erfcx, erfinv


    const Collection{T} = Union{AbstractArray{<:T}, AbstractRange{<:T}, NTuple{N,T}} where N
    include("Math.jl")

    using Colors: Color, RGBX, RGB, N0f8
    include("Interpolations.jl")
    include("ArrayStats.jl")

    using IndirectArrays: IndirectArray
    include("Images.jl")
    include("Colormaps.jl")

    import Base.display
    include("Display.jl") # Custom pretty-printing

    using DelimitedFiles
    include("Import.jl")
end

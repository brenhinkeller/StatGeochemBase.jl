module StatGeochemBase

    using NaNStatistics
    using VectorizedStatistics

    # AVX vectorziation tools
    using LoopVectorization: @turbo
    using VectorizationBase: Vec, verf
    using SpecialFunctions: erf, erfc, erfcx, erfinv
    import SpecialFunctions.erf
    erf(x::Vec) = verf(x)

    const Collection{T} = Union{DenseArray{<:T}, AbstractRange{<:T}, NTuple{N,T}} where N
    include("Math.jl")

    using Colors: Color, RGBX, RGB, N0f8
    include("Interpolations.jl")
    include("ArrayStats.jl")

    using IndirectArrays: IndirectArray
    include("Images.jl")
    include("Colormaps.jl")

    import Base.display
    include("Display.jl") # Custom pretty-printing

end

module StatGeochemBase

    using NaNStatistics

    # AVX vectorziation tools
    using LoopVectorization
    using VectorizationBase
    using VectorizationBase: Vec

    using SpecialFunctions
    import SpecialFunctions.erf
    erf(x::Vec) = VectorizationBase.verf(x)
    include("Math.jl")
    include("Geochronology.jl")

    using Colors: Color, RGBX, RGB, N0f8
    include("Interpolations.jl")
    include("ArrayStats.jl")

    using IndirectArrays: IndirectArray
    include("Images.jl")
    include("Colormaps.jl")

    import Base.display
    include("Display.jl") # Custom pretty-printing

end

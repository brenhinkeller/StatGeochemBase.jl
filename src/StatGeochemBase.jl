module StatGeochemBase

    using Reexport
    @reexport using NaNStatistics

    # AVX vectorziation tools
    using LoopVectorization
    using VectorizationBase
    using VectorizationBase: Vec

    using SpecialFunctions
    import SpecialFunctions.erf
    erf(x::Vec) = VectorizationBase.verf(x)

    # Other dependencies
    using Interpolations

    include("Math.jl")
    include("ArrayStats.jl")

end

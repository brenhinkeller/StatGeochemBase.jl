module StatGeochemBase

    using Reexport
    @reexport using NaNStatistics

    using LoopVectorization
    using Interpolations
    include("ArrayStats.jl")

end

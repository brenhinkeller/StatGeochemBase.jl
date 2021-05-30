module StatGeochemBase

    using Reexport
    @reexport using NaNStatistics

    using LoopVectorization
    include("ArrayStats.jl")

end

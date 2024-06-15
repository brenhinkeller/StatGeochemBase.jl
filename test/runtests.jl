using StatGeochemBase
using NaNStatistics
using Test

@testset "Math" begin include("testMath.jl") end
@testset "Images" begin include("testImages.jl") end
@testset "Import" begin include("testImport.jl") end
@testset "ArrayStats" begin include("testArrayStats.jl") end
@testset "Interpolations" begin include("testInterpolations.jl") end

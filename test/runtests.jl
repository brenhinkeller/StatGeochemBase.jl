using StatGeochemBase
using Test

@testset "Math" begin include("testMath.jl") end
@testset "Images" begin include("testImages.jl") end
@testset "ArrayStats" begin include("testArrayStats.jl") end
@testset "Geochronology" begin include("testGeochronology.jl") end

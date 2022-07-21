## --- test ArrayStats.jl

    # Sorting, counting, matching
    A = rand(1:100., 100); B = sort(A)
    @test A[1:count_unique!(A)] == unique(B)
    @test findclosest(3.3:5.3, 1:10) == 3:5
    @test findclosest(3.6, 1:10) == 4
    @test findclosestbelow(3.3:5.3, 1:10) == 3:5
    @test findclosestabove(3.3:5.3, 1:10) == 4:6
    @test findmatches(40:60, 1:100) == 40:60
    @test findmatches(50, 1:100) == 50
    @test findnth(fill(true, 50), 25) == 25

    # Interpolation
    @test cntr(0:2:100) == 1:2:99

    # Integration
    @test trapezoidalquadrature(1:10, fill(1,10)) == 9
    @test trapz(collect(1:10.), ones(10)) == 9
    @test midpointquadrature(1:10, ones(10)) == 10

    # Distributions
    A = draw_from_distribution(ones(100), 10000)::AbstractArray
    @test length(A) == 10000
    @test isapprox(nanmean(A), 0.5, atol=0.08)
    @test isapprox(nanstd(A), sqrt(1/12), atol=0.08)

    # Strings
    @test contains("JuliaLang is pretty cool!", "Julia")
    @test !contains("JuliaLang is pretty cool!", "julia")

    @test containsi("JuliaLang is pretty cool!", "Julia")
    @test containsi("JuliaLang is pretty cool!", "julia")
    @test !containsi("JuliaLang is pretty cool!", "tomatoes")


## ---

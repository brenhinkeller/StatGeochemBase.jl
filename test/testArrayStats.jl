## --- test ArrayStats.jl

    # Sorting, counting, matching
    A = rand(1:100.,100); B = sort(A)
    @test A[1:count_unique!(A)] == unique(B)
    @test findclosest(3.3:5.3,1:10) == 3:5
    @test findclosestbelow(3.3:5.3,1:10) == 3:5
    @test findclosestabove(3.3:5.3,1:10) == 4:6
    @test findmatches(40:60,1:100) == 40:60
    @test findnth(fill(true,50), 25) == 25

    # Interpolation
    @test linterp1(1:10,21:30,5:0.5:6) == [25.0, 25.5, 26.0]
    @test linterp1s(10:-1:1,21:30,5:0.5:6) == [26.0, 25.5, 25.0]
    @test linterp_at_index(1:100,10) == 10
    @test cntr(0:2:100) == 1:2:99

    # Integration
    @test trapezoidalquadrature(1:10,fill(1,10)) == 9
    @test trapz(collect(1:10.),ones(10)) == 9
    @test midpointquadrature(1:10,ones(10)) == 10

    # Distributions
    A = draw_from_distribution(ones(100), 10000)::AbstractArray
    @test length(A) == 10000
    @test isapprox(nanmean(A), 0.5, atol=0.08)
    @test isapprox(nanstd(A), sqrt(1/12), atol=0.08)

## ---

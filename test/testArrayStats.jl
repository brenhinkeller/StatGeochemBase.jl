## --- test ArrayStats.jl

    # Type wrangling
    a = Any[false, 0, 1.0]
    @test unionize(a) isa Vector{Union{Bool, Int, Float64}}
    @test unionize(a) == a
    @test unionize(1:10) === 1:10

    # Copying
    src = rand(100)
    t = src .< 0.5
    dest = fill(NaN, count(t))
    copyat!(dest, src, t)
    @test dest == src[t]

    reversecopyat!(dest, src, t)
    @test dest == reverse!(src[t])

    # Sorting, counting, matching
    A = rand(1:100., 100); B = sort(A)
    @test A[1:count_unique!(A)] == unique(B)
    @test findnth(fill(true, 50), 25) == 25
    @test findmatches(40:60, 1:100) == 40:60
    @test findmatches(50, 1:100) == 50

    @test findclosest(3.6, 1:10) == 4
    @test findclosest(-1, 1:10) == 1
    @test findclosest(11, 1:10) == 10
    @test findclosest(3.6, 10:-1:1) == 7
    @test findclosest(-1, 10:-1:1) == 10
    @test findclosest(11, 10:-1:1) == 1
    @test findclosest(3.6, [10, 3, 8, 6, 9, 2, 4, 7, 5, 1]) == 7
    @test findclosest(-1, [10, 3, 8, 6, 9, 2, 4, 7, 5, 1]) == 10
    @test findclosest(11, [10, 3, 8, 6, 9, 2, 4, 7, 5, 1]) == 1
    @test findclosest(3.3:5.3, 1:10) == 3:5
    @test findclosest(3.3:5.3, 10:-1:1) == 8:-1:6
    @test findclosest(3.3:5.3, [10, 3, 8, 6, 9, 2, 4, 7, 5, 1]) == [2,7,9]

    @test findclosestbelow(3.6, 1:10) == 3
    @test findclosestbelow(3.6, 10:-1:1) == 8
    @test findclosestbelow(3.6, (1, 3, 5, 10, 4, 7, 6, 8, 9)) == 2
    @test findclosestbelow(3.6, [1, 3, 5, 10, 4, 7, 6, 8, 9]) == 2
    @test findclosestbelow(3.6, [10, 3, 8, 6, 9, 2, 4, 7, 5, 1]) == 2
    @test findclosestbelow(-1, 1:10) == 0
    @test findclosestbelow(-1, 10:-1:1) == 0
    @test findclosestbelow(-1, [10, 3, 8, 6, 9, 2, 4, 7, 5, 1]) == 0
    @test findclosestbelow(11, 1:10) == 10
    @test findclosestbelow(11, 10:-1:1) == 1
    @test findclosestbelow(11, [10, 3, 8, 6, 9, 2, 4, 7, 5, 1]) == 1
    @test findclosestbelow(3.3:5.3, 1:10) == 3:5
    @test findclosestbelow((3.3:5.3...,), 1:10) == 3:5
    @test findclosestbelow(3.3:5.3, 10:-1:1) == 11 .- (3:5)
    @test findclosestbelow((3.3:5.3...,), 10:-1:1) == 11 .- (3:5)

    @test findclosestabove(3.6, 1:10) == 4
    @test findclosestabove(3.6, 10:-1:1) == 7
    @test findclosestabove(3.6, (1, 3, 5, 10, 4, 7, 6, 8, 9)) == 5
    @test findclosestabove(3.6, [1, 3, 5, 10, 4, 7, 6, 8, 9]) == 5
    @test findclosestabove(3.6, [10, 3, 8, 6, 9, 2, 4, 7, 5, 1]) == 7
    @test findclosestabove(11, 1:10) == 11
    @test findclosestabove(11, 10:-1:1) == 11
    @test findclosestabove(11, [10, 3, 8, 6, 9, 2, 4, 7, 5, 1]) == 11
    @test findclosestabove(0, 1:10) == 1
    @test findclosestabove(0, 10:-1:1) == 10
    @test findclosestabove(0, [10, 3, 8, 6, 9, 2, 4, 7, 5, 1]) == 10
    @test findclosestabove(3.3:5.3, 1:10) == 4:6
    @test findclosestabove((3.3:5.3...,), 1:10) == 4:6
    @test findclosestabove(3.3:5.3, 10:-1:1) == 11 .- (4:6)
    @test findclosestabove((3.3:5.3...,), 10:-1:1) == 11 .- (4:6)

    x = fill(1, 50)
    @test findclosestunequal(x, 25) == 25
    x[end] = 2
    @test findclosestunequal(x, 25) == 50
    x[1] = 0
    @test findclosestunequal(x, 25) == 1

    # Interpolation
    @test cntr(0:2:100) == 1:2:99
    @test stepify([1,3,2]) == [1, 1, 3, 3, 2, 2]
    @test stepifyedges([1,2,3,4]) == [1, 2, 2, 3, 3, 4]

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

    # MCMC inversion on a precalculated surface
    xaxis = yaxis = -5:0.01:5
    llmatrix = [-((x+1)^2/0.5^2 + (y-1)^2/0.5^2) for y in yaxis, x in xaxis];
    xdist, ydist, lldist, acceptancedist = mcmc_surface(xaxis,yaxis,llmatrix);

    @test nanmean(xdist) ≈ -1 atol=0.1
    @test nanmean(ydist) ≈ 1 atol=0.1
    @test nanstd(xdist) ≈ 0.5 atol=0.15
    @test nanstd(ydist) ≈ 0.5 atol=0.15
    @test nanmean(acceptancedist) ≈ 0.32 atol=0.2

## ---

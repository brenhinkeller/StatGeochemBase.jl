## --- test Interpolations.jl

    # Interpolation
    @test linterp1(1:10, 1:10, 5.5) == 5.5
    @test linterp1(1:10, collect(1:10.), 3:7) == 3:7
    @test linterp1(1:10,21:30,5:0.5:6) == [25.0, 25.5, 26.0]
    @test linterp1s(10:-1:1,21:30,5:0.5:6) == [26.0, 25.5, 25.0]
    @test linterp_at_index(1:100,10) == 10

    # Extrapolation
    @test linterp1(1:10, 1:10, 15) == 15 # Default is to extrapolate
    @test linterp1(1:10, 1:10, 15, extrapolate=-5) == -5
    @test linterp1(1:10, 1:10, 5, extrapolate=-5) == 5
    @test isnan(linterp1(1:10, 1:10, 15, extrapolate=NaN))
    @test linterp1(1:10,1:10,0:11) == 0:11 # Default is to extrapolate
    @test linterp1(1:10,1:10,0:11, extrapolate=:Linear) == 0:11
    @test linterp1(1:10,1:10,0.5:10.5, extrapolate=:Linear) == 0.5:10.5
    @test linterp1(1:10,1:10,0.5:10.5, extrapolate=-5) == [-5; 1.5:9.5; -5]
    @test all(linterp1(1:10,1:10,0.5:10.5, extrapolate=NaN) .=== [NaN; 1.5:9.5; NaN])
    @test isnan(linterp_at_index(1:100,-10))
    @test linterp_at_index(1:100,-10, 0) == 0

    # In-place
    xq = 3:7
    @test linterp1!(similar(xq, Float64), 1:10, collect(1:10.), xq) == 3:7
    xq = 5:0.5:6
    @test linterp1!(similar(xq), 1:10, 21:30, xq) == [25.0, 25.5, 26.0]
    xq = 5:0.5:6

    @test linterp1s!(similar(xq), collect(1:10), collect(21:30), xq) == [25.0, 25.5, 26.0]
    @test linterp1s!(similar(xq), collect(10:-1:1), collect(21:30), xq) == [26.0, 25.5, 25.0]
    xq = 5:0.5:7
    @test linterp1s!(similar(xq), rand(Int,length(xq)), collect(10:-1:1), collect(21:30), xq) == [26.0, 25.5, 25.0, 24.5, 24]
    xq = 0:0.01:1
    x = rand(200)
    y = rand(200)
    yq = linterp1s(x, y, xq)
    @test linterp1s!(similar(xq), x, y, xq) ≈ yq

    xq = 0:11
    @test linterp1!(similar(xq, Float64), 1:10, 1:10, xq, extrapolate=:Linear) == 0:11
    xq = 0.5:10.5
    @test isequal(linterp1!(similar(xq), 1:10, 1:10, xq, extrapolate=NaN), [NaN; 1.5:9.5; NaN])

    # Test consistency of sorting against Base
    x = rand(10)*10
    y = rand(10)*10
    perm = sortperm(x)
    xs = x[perm]
    yx = y[perm]

    xq = 0:0.01:10
    yq = similar(xq)
    knot_index = rand(Int,length(xq))
    linterp1s!(yq, knot_index, x, y, xq)
    @test yq == linterp1(xs, yx, xq)

    x = rand(1000)*10
    y = rand(1000)*10
    perm = sortperm(x)
    xs = x[perm]
    yx = y[perm]

    xq = 0:0.01:10
    yq = similar(xq)
    knot_index = rand(Int,length(xq))
    linterp1s!(yq, knot_index, x, y, xq)
    @test yq == linterp1(xs, yx, xq)

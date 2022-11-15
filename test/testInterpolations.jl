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
    @test linterp1s!(similar(xq), 1:10, 21:30, xq) == [25.0, 25.5, 26.0]
    @test linterp1s!(similar(xq), collect(10:-1:1), collect(21:30), xq) == [26.0, 25.5, 25.0]
    xq = 5:0.5:7
    @test linterp1s!(similar(xq), collect(10:-1:1), collect(21:30), xq; knot_index=ones(Int,5), sort_index=ones(Int,10)) == [26.0, 25.5, 25.0, 24.5, 24]
    xq = 0:11
    @test linterp1!(similar(xq, Float64), 1:10, 1:10, xq, extrapolate=:Linear) == 0:11
    xq = 0.5:10.5
    @test isequal(linterp1!(similar(xq), 1:10, 1:10, xq, extrapolate=NaN), [NaN; 1.5:9.5; NaN])

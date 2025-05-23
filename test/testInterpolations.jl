## --- Test 1-D interpolations

    # Interpolation
    @test linterp1(1:10, 1:10, 1, extrapolate=NaN) == 1
    @test linterp1(1:10, 1:10, 10, extrapolate=NaN) == 10
    @test linterp1(1:10, 1:10, 5.5) == 5.5
    @test linterp1(1:10, 1:10, 1:10, extrapolate=NaN) == 1:10
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

## --- Test 2-D interpolations

    x = 1:3
    y = 1:4
    z = y*x'
    @test linterp2(x,y,z,1.5,1.5) ≈ 2.25
    @test linterp2(x,y,z,1.5,1.5, extrapolate=NaN) ≈ 2.25
    @test linterp2(x,y,z,1.5,1.5, extrapolate=:Bilinear) ≈ 2.25
    @test linterp2(x,y,z,2.5,3.5) ≈ 8.75
    @test linterp2(x,y,z,1,-10,extrapolate=:Bilinear) ≈ -10.0
    @test linterp2(x,y,z,-10,-10,extrapolate=:Bilinear) ≈ 100
    @test isnan(linterp2(x,y,z,1,-10,extrapolate=NaN))
    @test isnan(linterp2(x,y,z,-10,-10,extrapolate=NaN))

    z = rand(4,3)
    @test linterp2(x,y,z,1.5,1.5) ≈ linterp2(x,y,z,1.5,1.5, extrapolate=NaN) ≈ nanmean(z[1:2,1:2])
    @test linterp2(x,y,z,2.5,3.5) ≈ linterp2(x,y,z,2.5,3.5, extrapolate=NaN) ≈ nanmean(z[3:4,2:3])
    
    xq = rand()*2+1
    yq = rand()*3+1
    @test linterp2(x,y,z,xq,yq) ≈ linterp2(x,y,z,xq,yq, extrapolate=NaN)
    @test linterp2(x,y,z,xq,yq) ≈ linterp2(x,y,z,xq,yq, extrapolate=NaN)

## --- End of File
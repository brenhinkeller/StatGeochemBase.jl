## --- Test interpolation-related utilities
    v1 = Base.OneTo(100)
    v2 = 1:100
    v3 = 1:100.
    v4 = collect(v3)

    x = 110*rand(1000) .- 5
    ix = zeros(Int, length(x))
    @test searchsortedfirst.(Ref(v1),x) == StatGeochemBase.searchsortedfirst_vec!(copy(ix), v1, x) == StatGeochemBase.searchsortedfirst_vec!(copy(ix), v2, x) == StatGeochemBase.searchsortedfirst_vec!(copy(ix), v3, x) == StatGeochemBase.searchsortedfirst_vec!(copy(ix), v4, x) 

    v2 = 10:90
    v3 = 10:90.
    v4 = collect(v3)
    @test searchsortedfirst.(Ref(v2),x) == StatGeochemBase.searchsortedfirst_vec!(copy(ix), v2, x) == StatGeochemBase.searchsortedfirst_vec!(copy(ix), v3, x) == StatGeochemBase.searchsortedfirst_vec!(copy(ix), v4, x) 

## --- Test 1-D interpolations

    # Interpolation
    @test linterp1(1:10, 1:10, 1, extrapolate=NaN) == 1
    @test linterp1(1:10, 1:10, 10, extrapolate=NaN) == 10
    @test linterp1(1:10, 1:10, 5.5) == 5.5
    @test linterp1(1:10, 1:10, 1:10, extrapolate=NaN) == 1:10
    @test linterp1(1:10, collect(1:10.), 3:7) == 3:7
    @test linterp1(1:10,21:30,5:0.5:6) == [25.0, 25.5, 26.0]
    @test linterp1s(10:-1:1,21:30,5:0.5:6) == [26.0, 25.5, 25.0]
    @test linterp1(1:10, 1:10, 1, interpolate=:nearest, extrapolate=NaN) == 1
    @test linterp1(1:10, 1:10, 10, interpolate=:nearest, extrapolate=NaN) == 10
    @test linterp1(1:10, 1:10, 5.4, interpolate=:nearest) == 5
    @test linterp1(1:10, 1:10, 1:10, interpolate=:nearest, extrapolate=NaN) == 1:10
    @test linterp1(1:10, collect(1:10.), 3:7, interpolate=:nearest,) == 3:7
    @test linterp1(1:10,21:30,5:0.5:6, interpolate=:nearest,) == [25., 25., 26.]
    @test linterp1s(10:-1:1,21:30,5:0.5:6, interpolate=:nearest,) == [26., 26., 25.]
    @test linterp_at_index(1:100,10) == 10

    # Extrapolation
    @test linterp1(1:10, 1:10, 15) == 15 # Default is to extrapolate
    @test linterp1(1:10, 1:10, 15, extrapolate=-5) == -5
    @test linterp1(1:10, 1:10, 5.5, extrapolate=-5) == 5.5
    @test linterp1(1:10, 1:10, 15, extrapolate=:nearest) == 10
    @test linterp1(1:10, 1:10, 5.5, extrapolate=:nearest) == 5.5
    @test isnan(linterp1(1:10, 1:10, 15, extrapolate=NaN))
    @test linterp1(1:10, 1:10, 15, interpolate=:nearest, extrapolate=-5) == -5
    @test linterp1(1:10, 1:10, 15, interpolate=:nearest, extrapolate=:nearest) == 10
    @test linterp1(1:10, 1:10, 5.5, interpolate=:nearest, extrapolate=:nearest) == 5
    @test isnan(linterp1(1:10, 1:10, 15, extrapolate=NaN))
    @test linterp1(1:10,1:10,0:11) == 0:11 # Default is to extrapolate
    @test linterp1(1:10,1:10,0:11, extrapolate=:linear) == 0:11
    @test linterp1(1:10,1:10,0:11, extrapolate=:nearest) == [1; 1:10; 10]
    @test linterp1(1:10,1:10,0.5:10.5, extrapolate=:linear) == 0.5:10.5
    @test linterp1(1:10,1:10,0.5:10.5, extrapolate=:nearest) == [1; 1.5:9.5; 10]
    @test linterp1(1:10,1:10,0.5:10.5, extrapolate=-5) == [-5; 1.5:9.5; -5]
    @test linterp1(1:10,1:10,0.5:10.5, extrapolate=NaN) ≈ [NaN; 1.5:9.5; NaN] nans=true
    @test linterp1(1:10,1:10,0:11, interpolate=:nearest, extrapolate=:linear) == 0:11
    @test linterp1(1:10,1:10,0:11, interpolate=:nearest, extrapolate=:nearest) == [1; 1:10; 10]
    @test linterp1(1:10,1:10,0.5:10.5, interpolate=:nearest, extrapolate=:linear) == [0.5; 1:9; 10.5]
    @test linterp1(1:10,1:10,0.5:10.5, interpolate=:nearest, extrapolate=:nearest) == [1; 1:9; 10]
    @test linterp1(1:10,1:10,0.5:10.5, interpolate=:nearest, extrapolate=-5) == [-5; 1:9; -5]
    @test linterp1(1:10,1:10,0.5:10.5, interpolate=:nearest, extrapolate=NaN) ≈ [NaN; 1:9; NaN] nans=true
    @test isnan(linterp_at_index(1:100,-10))
    @test linterp_at_index(1:100,-10, 0) == 0

    # In-place
    xq = 3:7
    @test linterp1!(similar(xq, Float64), 1:10, collect(1:10.), xq) == 3:7
    @test linterp1!(similar(xq, Float64), 1:10, collect(1:10.), xq, interpolate=:nearest) == 3:7

    xq = 5:0.5:6
    @test linterp1!(similar(xq), 1:10, 21:30, xq) == [25.0, 25.5, 26.0]
    @test linterp1!(similar(xq), 1:10, 21:30, xq, interpolate=:nearest) == [25., 25., 26.]

    xq = 5:0.5:6
    @test linterp1s!(similar(xq), collect(1:10), collect(21:30), xq) == [25.0, 25.5, 26.0]
    @test linterp1s!(similar(xq), collect(10:-1:1), collect(21:30), xq) == [26.0, 25.5, 25.0]
    @test linterp1s!(similar(xq), collect(1:10), collect(21:30), xq, interpolate=:nearest) == [25., 25., 26.]
    @test linterp1s!(similar(xq), collect(10:-1:1), collect(21:30), xq, interpolate=:nearest) == [26., 26., 25.]

    xq = 5:0.5:7
    @test linterp1s!(similar(xq), rand(Int,length(xq)), collect(10:-1:1), collect(21:30), xq) == [26.0, 25.5, 25.0, 24.5, 24]
    @test linterp1s!(similar(xq), rand(Int,length(xq)), collect(10:-1:1), collect(21:30), xq, interpolate=:nearest) == [26., 26., 25., 25., 24]

    xq = 0:0.01:1
    x = rand(200)
    y = rand(200)
    yq = linterp1s(x, y, xq, interpolate=:linear)
    @test linterp1s!(similar(xq), x, y, xq, interpolate=:linear) ≈ yq
    yq = linterp1s(x, y, xq, interpolate=:nearest)
    @test linterp1s!(similar(xq), x, y, xq, interpolate=:nearest) ≈ yq

    xq = 0:11
    @test linterp1!(similar(xq, Float64), 1:10, 1:10, xq, interpolate=:linear) == 0:11
    @test linterp1!(similar(xq, Float64), 1:10, 1:10, xq, interpolate=:nearest) == [1; 1:10; 10]

    xq = 0.5:10.5
    @test isequal(linterp1!(similar(xq), 1:10, 1:10, xq, interpolate=:linear, extrapolate=NaN), [NaN; 1.5:9.5; NaN])
    @test isequal(linterp1!(similar(xq), 1:10, 1:10, xq, interpolate=:nearest, extrapolate=NaN), [NaN; 1:9; NaN])

    # Test consistency of sorting against Base
    x = rand(50)*10
    y = rand(50)*10
    perm = sortperm(x)
    xs = x[perm]
    yx = y[perm]

    xq = 0:0.01:10
    yq = similar(xq)
    linterp1s!(yq, x, y, xq)
    @test yq ≈ linterp1(xs, yx, xq)
    linterp1s!(yq, x, y, xq, extrapolate=NaN)
    @test yq ≈ linterp1(xs, yx, xq, extrapolate=NaN) nans=true
    linterp1s!(yq, x, y, xq, interpolate=:nearest)
    @test yq ≈ linterp1(xs, yx, xq, interpolate=:nearest)
    linterp1s!(yq, x, y, xq, interpolate=:nearest, extrapolate=NaN)
    @test yq ≈ linterp1(xs, yx, xq, interpolate=:nearest, extrapolate=NaN) nans=true

    yqk = similar(xq)
    knot_index = rand(Int,length(xq))
    linterp1s!(yqk, knot_index, x, y, xq)
    @test yqk ≈ linterp1(xs, yx, xq)
    linterp1s!(yqk, knot_index, x, y, xq, extrapolate=NaN)
    @test yqk ≈ linterp1(xs, yx, xq, extrapolate=NaN) nans=true
    linterp1s!(yqk, knot_index, x, y, xq, interpolate=:nearest)
    @test yqk ≈ linterp1(xs, yx, xq, interpolate=:nearest)
    linterp1s!(yqk, knot_index, x, y, xq, interpolate=:nearest, extrapolate=NaN)
    @test yqk ≈ linterp1(xs, yx, xq, interpolate=:nearest, extrapolate=NaN) nans=true

    x = rand(5000)*10
    y = rand(5000)*10
    perm = sortperm(x)
    xs = x[perm]
    yx = y[perm]

    xq = 0:0.01:10
    yq = similar(xq)
    linterp1s!(yq, x, y, xq)
    @test yq ≈ linterp1(xs, yx, xq)
    linterp1s!(yq, x, y, xq, interpolate=:nearest)
    @test yq ≈ linterp1(xs, yx, xq, interpolate=:nearest)

    yqk = similar(xq)
    knot_index = rand(Int,length(xq))
    linterp1s!(yqk, knot_index, x, y, xq)
    @test yqk ≈ linterp1(xs, yx, xq)
    linterp1s!(yqk, knot_index, x, y, xq, interpolate=:nearest)
    @test yqk ≈ linterp1(xs, yx, xq, interpolate=:nearest)

## --- Test 2-D interpolations

    x = 1:3
    y = 1:4
    z = y*x'
    @test linterp2(x,y,z,1.5,1.5) ≈ 2.25
    @test linterp2(x,y,z,1.5,1.5, extrapolate=NaN) ≈ 2.25
    @test linterp2(x,y,z,1.5,1.5, interpolate=:bilinear) ≈ 2.25
    @test linterp2(x,y,z,1.6,1.6, interpolate=:nearest) ≈ 4
    @test linterp2(x,y,z,2.5,3.5) ≈ 8.75
    @test linterp2(x,y,z,1,-10,extrapolate=:bilinear) ≈ -10.0
    @test linterp2(x,y,z,-10,-10,extrapolate=:bilinear) ≈ 100
    @test linterp2(x,y,z,10,10,extrapolate=:bilinear) ≈ 100
    @test linterp2(x,y,z,1,-10,extrapolate=:nearest) ≈ 1
    @test linterp2(x,y,z,-10,-10,extrapolate=:nearest) ≈ 1
    @test linterp2(x,y,z,10,10,extrapolate=:nearest) ≈ 12
    @test isnan(linterp2(x,y,z,1,-10,extrapolate=NaN))
    @test isnan(linterp2(x,y,z,-10,-10,extrapolate=NaN))

    z = rand(4,3)
    @test linterp2(x,y,z,1.5,1.5) ≈ linterp2(x,y,z,1.5,1.5, extrapolate=NaN) ≈ nanmean(z[1:2,1:2])
    @test linterp2(x,y,z,2.5,3.5) ≈ linterp2(x,y,z,2.5,3.5, extrapolate=NaN) ≈ nanmean(z[3:4,2:3])
    
    xq = rand()*2+1
    yq = rand()*3+1
    @test linterp2(x,y,z,xq,yq) ≈ linterp2(x,y,z,xq,yq, extrapolate=NaN)
    @test linterp2(x,y,z,xq,yq) ≈ linterp2(x,y,z,xq,yq, extrapolate=NaN)


    # Test 2d interpolation with a matrix of xq and yq values
    x = 1:3
    y = 1:4
    z = y*x'
    xq = repeat(x', 4, 1)
    yq = repeat(y, 1, 3)
    @test linterp2(x,y,z,xq,yq) ≈ z
    @test linterp2(x,y,z,xq,yq, interpolate=:nearest) ≈ z

## --- End of File
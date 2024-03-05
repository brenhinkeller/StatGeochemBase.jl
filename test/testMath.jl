## --- Special functions

    @test isapprox(fast_inv_sqrt(5.0), 1/sqrt(5.0), atol=1e-6)
    @test isapprox(fast_inv_sqrt(5f0), 1/sqrt(5f0), atol=1e-6)
    @test nearest(Int64, 3.3) === 3
    @test nearest(Float64, 1//3) === 1/3

## --- Significant figures

    @test sigdigits(1/3) === 16
    @test sigdigits(0.11) === 2
    @test sigdigits(0.111) === 3
    @test sigdigits(0.1111) === 4
    @test sigdigits(0.11111) === 5
    @test sigdigits(0.111111) === 6
    @test sigdigits(0.1111111) === 7
    @test sigdigits(0.11111111) === 8
    @test sigdigits(0.111111111) === 9
    @test sigdigits(0.1111111111) === 10
    @test sigdigits(0.11111111111) === 11
    @test sigdigits(0.111111111111) === 12
    @test sigdigits(0.1111111111111) === 13
    @test sigdigits(0.11111111111111) === 14
    @test sigdigits(0.111111111111111) === 15
    @test sigdigits(0.1111111111111111) === 16

    for T in (BigInt, Int64, Int32, Int16, Int8, UInt64, UInt32, UInt16, UInt8)
        @test sigdigits(T(100)) === 1
        @test sigdigits(T(101)) === 3
        @test leastsigfig(T(100)) === 100.
        @test leastsigfig(T(101)) === 1.
    end
    @test sigdigits(big"1000.") === sigdigits(1000.) === sigdigits(Float32(1000.)) === sigdigits(Float16(1000.)) === 1
    @test sigdigits(big"1000.5") === sigdigits(1000.5) === sigdigits(Float32(1000.5)) === 5
    @test leastsigfig(big"1000.") === leastsigfig(1000.) === leastsigfig(Float32(1000.)) === leastsigfig(Float16(1000.)) === 1000.
    @test leastsigfig(big"1000.5") === leastsigfig(1000.5) === leastsigfig(Float32(1000.5)) === 0.1

## --- Linear regression

    @test linreg(1:10, 1:10) ≈ [0, 1]

## --- Distributions

    @test normpdf.(0, 1,-1:1) ≈ [0.24197072451914337, 0.3989422804014327, 0.24197072451914337]
    @test normpdf.(1:10, 1:10, 1:10) ≈ normpdf(collect.((1:10, 1:10, 1:10))...)

    @test normpdf_ll.(0,1,-5:5) == -(-5:5).^2/2
    r = collect(-5:5)
    @test normpdf_ll(0,1,r) == normpdf_ll(0,ones(11),r) == normpdf_ll(zeros(11),ones(11),r) == sum(normpdf_ll.(0,1,r))
    @test normpdf_ll(ones(10),1,collect(1:10)) == normpdf_ll(collect(1:10),1,ones(10)) ≈ -142.5 # Test for symmetry

    @test normcdf(1,1,1) == 0.5
    result = zeros(5)
    normcdf!(result, 0, 1, -2:2)
    @test result ≈ normcdf(0,1,-2:2) ≈ normcdf.(0,1,-2:2) ≈ [0.02275013194817921, 0.15865525393145707, 0.5, 0.8413447460685429, 0.9772498680518208]
    @test normcdf.(1:10, 1:10, 1:10) == normcdf(collect.((1:10, 1:10, 1:10))...) == fill(0.5, 10)

    @test normcdf_ll.(0,1,-5:5) ≈ [-15.064998393988725, -10.360101486527292, -6.607726221510349, -3.7831843336820317, -1.841021645009264, -0.6931471805599453, -0.17275377902344985, -0.023012909328963486, -0.0013508099647481923, -3.1671743377489226e-5, -2.866516129637633e-7]
    r = collect(-5:5)
    @test normcdf_ll(r) ≈ -38.54732871798976
    @test normcdf_ll(r) == normcdf_ll(0,1,r) == normcdf_ll(0,ones(11),r) == normcdf_ll(zeros(11),ones(11),r) == sum(normcdf_ll.(0,1,r))
    @test normcdf_ll(zeros(10),1,collect(1:10)) ==  normcdf_ll(-collect(1:10),1,zeros(10)) ≈ -0.19714945770002004 # Test for symmetry
    @test normcdf_ll!(float.(r)) ≈ -38.54732871798976
    @test normcdf_ll!(float.(r)) == normcdf_ll!(0,1,float.(r)) == normcdf_ll!(0,ones(11),float.(r)) == normcdf_ll!(zeros(11),ones(11),float.(r))
    @test normcdf_ll!(zeros(10),1,collect(1:10.)) == normcdf_ll!(-collect(1:10),1,zeros(10)) ≈ -0.19714945770002004 # Test for symmetry

    @test normproduct(0,1,0,1) === normpdf(0,sqrt(2),0) === 0.28209479177387814
    @test normproduct_ll(0,1,0,1) === normpdf_ll(0,1,0) === 0.0

    @test [-2,0,2] ≈ norm_quantile.([0.022750131948, 0.5, 0.977249868052])
    @test norm_quantile.(0:0.25:1) ≈ [-Inf, -0.6744897501960818, 0.0, 0.6744897501960818, Inf]

    @test isapprox(norm_width(390682215445)/2, 7, atol=1e-5)


## -- Geometry

    @test inpolygon([-1,0,1,0],[0,1,0,-1],[0,0])
    @test !inpolygon([-1,0,1,0],[0,1,0,-1],[0,10])
    @test inpolygon([-1,0,1,0],[0,1,0,-1],prevfloat.([0.5,0.5]))
    @test !inpolygon([-1,0,1,0],[0,1,0,-1],nextfloat.([0.5,0.5]))
    @test inpolygon([-1,1,1,-1],[1,1,-1,-1],(0,0))
    @test !inpolygon([-1,1,1,-1],[1,1,-1,-1],(1.1,1))

    i,j = find_grid_inpolygon(-1.5:1/3:1.5, -1.5:1/3:1.5, [-.75,.75,.75,-.75],[.75,.75,-.75,-.75])
    @test sort([i j], dims=2) == [4 4; 4 5; 4 6; 4 7; 4 5; 5 5; 5 6; 5 7; 4 6; 5 6; 6 6; 6 7; 4 7; 5 7; 6 7; 7 7]

    @test arcdistance(0,100,[30,0,0],[100,100,95]) ≈ [30,0,5]

    @test minarcdistance(0,100,[30,0,0],[100,100,95]) ≈ fill(0)
    @test minarcdistance([1,0,1,2,3,4],[101,100,100,100,100,100],[30,0,0],[100,100,95]) ≈ [1.414177660951948,0,1,2,3,4]


## ---

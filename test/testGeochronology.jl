## --- Weighted means

    # Plain vanilla weighted mean
    @test awmean([0,1,2,3],[1,1,1,1]) == (1.5, 0.5, 5/3)
    @test awmean(1:10,ones(10)) == (5.5, 0.31622776601683794, 9.166666666666666)

    # Weighted mean with geochronogists' MSWD-correction to uncertainty
    @test gwmean([0,1,2,3],[1,1,1,1]) == (1.5, 0.6454972243679028, 5/3)
    @test gwmean(1:10,ones(10)) == (5.5, 0.9574271077563381, 9.166666666666666)

    # Mean Square of Weighted Deviates (aka reduced chi-squared)
    @test MSWD(0:10, ones(11)) == 11
    @test MSWD([0,1,2],[1,1,1]) == 1.0

## --- York fit

    x = [0.9304, 2.2969, 2.8047, 3.7933, 5.3853, 6.1995, 6.7479, 8.1856, 8.7423, 10.2588]
    y = [0.8742, 2.1626, 3.042, 3.829, 5.0116, 5.5614, 6.7675, 7.8856, 9.6414, 10.4955]
    fobj = yorkfit(x, ones(10)/4, y, ones(10)/4)
    @test fobj.intercept ≈  -0.23498964673701916
    @test fobj.intercept_sigma ≈ 0.02250863813481163
    @test fobj.slope ≈ 1.041124018512526
    @test fobj.slope_sigma ≈ 0.0035683808205783673
    @test fobj.mswd ≈ 1.1419901440278089

## ---

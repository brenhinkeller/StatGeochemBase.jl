## --- Images.jl

   using Random
   rng = Xoshiro(1234)
   
   nsims = 1000
   xpaths = rand(rng, 10, nsims)
   ypaths = rand(rng, 10, nsims)
   xpaths[1,:] .= 0
   xpaths[end,:] .= 1
   imgcounts, xq, yq = image_from_paths(xpaths, ypaths; xresolution=100, yresolution=50, xrange=(0,1), yrange=(0,1), method=:interpolate)
   @test sum(imgcounts, dims=1) == fill(1000, 1, 100)
   @test sum(imgcounts, dims=2) == [205; 674; 745; 1029; 1219; 1383; 1580; 1890; 1935; 1993; 2146; 2068; 2241; 2388; 2278; 2384; 2505; 2534; 2618; 2684; 2681; 2766; 2758; 2860; 2734; 2781; 2827; 2851; 2758; 2730; 2797; 2657; 2591; 2538; 2478; 2265; 2284; 2210; 2132; 2032; 1920; 1782; 1721; 1617; 1267; 1082; 968; 705; 515; 194;;]
   @test xq == 0.005:0.01:0.995
   @test yq == 0.01:0.02:0.99

   nsims = 1500
   xpaths = rand(rng, 10, nsims)
   ypaths = rand(rng, 10, nsims)
   xpaths[1,:] .= 0
   xpaths[end,:] .= 1
   imgcounts, xq, yq = image_from_paths(xpaths, ypaths; xresolution=100, yresolution=50, xrange=(0,1), yrange=(0,1), method=:interpolate)
   @test xq == 0.005:0.01:0.995
   @test yq == 0.01:0.02:0.99
   @test sum(imgcounts, dims=1) == fill(1500, 1, 100)
   @test sum(imgcounts, dims=2) == [289; 898; 1211; 1584; 1784; 2132; 2383; 2532; 2643; 2835; 3066; 3048; 3284; 3503; 3575; 3654; 3770; 3795; 3865; 3864; 3861; 3919; 4162; 4207; 4111; 4053; 4038; 4264; 4157; 4188; 4095; 4127; 3977; 4083; 4004; 3807; 3672; 3530; 3336; 3034; 2907; 2733; 2465; 2201; 1967; 1655; 1499; 1157; 767; 309;;]

   imgcounts = zeros(12,24)
   image_from_paths!(imgcounts, xpaths, ypaths; xrange=(0,1), yrange=(0,1), method=:interpolate)
   @test sum(imgcounts, dims=1) == fill(1500, 1, 24)
   @test sum(imgcounts, dims=2) == [1009.0; 2215.0; 3026.0; 3520.0; 3862.0; 4127.0; 4061.0; 4138.0; 3863.0; 3104.0; 2138.0; 937.0;;]

   imgcounts = zeros(12,24)
   image_from_paths!(imgcounts, xpaths, ypaths; xrange=(0,1), yrange=(0,1), method=:nearest)
   @test sum(imgcounts, dims=1) == [4128.0 2842.0 2741.0 2704.0 2571.0 2486.0 2503.0 2481.0 2530.0 2573.0 2634.0 2593.0 2605.0 2589.0 2532.0 2417.0 2610.0 2539.0 2430.0 2480.0 2592.0 2506.0 2862.0 4359.0]
   # println(sum(imgcounts, dims=1))
   @test sum(imgcounts, dims=2) == [1566.0; 3774.0; 5320.0; 6612.0; 7290.0; 7763.0; 7778.0; 7640.0; 6855.0; 5490.0; 3661.0; 1558.0;;]
   # println(sum(imgcounts, dims=2))

   using Colors: Color

   cmap = resize_colormap(viridis, 10)
   @test length(cmap) == 10
   @test isa(cmap, Array{<:Color,1})

   matrix = rand(10,10)
   # Specifiying limits
   img1 = imsc(matrix, viridis, 0, 1)
   @test isa(img1, Array{<:Color,2})
   img2 = imsci(matrix, viridis, 0, 1)
   @test isa(img2, AbstractArray{<:Color,2})
   @test all(img1 .== img2)

   # Auto-ranging
   img1 = imsc(matrix, viridis)
   @test isa(img1, Array{<:Color,2})
   img2 = imsci(matrix, viridis)
   @test isa(img2, AbstractArray{<:Color,2})
   @test all(img1 .== img2)

   # Other
   #@test display(colormaps) != NaN

## ---

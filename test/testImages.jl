## --- Images.jl

   using Random
   rng = Xoshiro(1234)
   
   nsims = 1000
   xpaths = rand(rng, 10, nsims)
   ypaths = rand(rng, 10, nsims)
   imgcounts, xq, yq = image_from_paths(xpaths, ypaths; xresolution=50, yresolution=50, xrange=(0,1), yrange=(0,1))

   @test sum(imgcounts, dims=1) == [715 999 1107 1268 1424 1502 1624 1636 1810 1837 1858 1925 1970 2041 2112 2059 2114 2139 2175 2175 2193 2200 2237 2265 2217 2169 2242 2303 2257 2225 2236 2187 2127 2130 2121 2060 2023 2015 1948 1870 1858 1755 1731 1643 1476 1410 1250 1103 916 713]
   @test sum(imgcounts, dims=2) == [715; 999; 1107; 1268; 1424; 1502; 1624; 1636; 1810; 1837; 1858; 1925; 1970; 2041; 2112; 2059; 2114; 2139; 2175; 2175; 2193; 2200; 2237; 2265; 2217; 2169; 2242; 2303; 2257; 2225; 2236; 2187; 2127; 2130; 2121; 2060; 2023; 2015; 1948; 1870; 1858; 1755; 1731; 1643; 1476; 1410; 1250; 1103; 916; 713;;]
   @test xq == 0.01:0.02:0.99
   @test yq == 0.01:0.02:0.99

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

## --- Custom display functions

    # Custom pretty-printing for colormaps
    function display(x::AllColormaps)
        println("AllColormaps:")
        for name in fieldnames(AllColormaps)
            println("  $name")
            display(getfield(x, name))
        end
    end

    # Custom pretty-printing for York fit results
    function display(x::YorkFit)
        print("York Fit y = a + bx:\n  intercept a: $(x.intercept) ± $(x.intercept_sigma) (1σ)")
        print("\n  slope b    : $(x.slope) ± $(x.slope_sigma) (1σ)\n  MSWD       : $(x.mswd)\n")
    end

## --- 

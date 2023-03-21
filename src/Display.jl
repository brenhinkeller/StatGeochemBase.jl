## --- Custom display functions

    # Custom pretty-printing for colormaps
    function display(x::AllColormaps)
        println("AllColormaps:")
        for name in fieldnames(AllColormaps)
            println("  $name")
            display(getfield(x, name))
        end
    end


## ---

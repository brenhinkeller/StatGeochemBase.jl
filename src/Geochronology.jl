## --- Weighted mean of an array

    """
    ```julia
    (wx, wσ, mswd) = awmean(x, σ)
    ```
    Weighted mean, absent the geochonologist's MSWD correction to uncertainty.
    """
    function awmean(x, σ)
        n = length(x)

        sum_of_values = sum_of_weights = χ2 = 0.0
        @inbounds @simd for i=1:n
            sum_of_values += x[i] / (σ[i]*σ[i])
            sum_of_weights += 1 / (σ[i]*σ[i])
        end
        wx = sum_of_values / sum_of_weights

        @inbounds @simd for i=1:n
            χ2 += (x[i] - wx) * (x[i] - wx) / (σ[i] * σ[i])
        end
        mswd = χ2 / (n-1)
        wσ = sqrt(1.0 / sum_of_weights)
        return wx, wσ, mswd
    end
    # function awmean(x::Array{<:Number}, σ::Array{<:Number})
    #     n = length(x)
    #
    #     sum_of_values = sum_of_weights = χ2 = 0.0
    #     @avx for i=1:n
    #         sum_of_values += x[i] / (σ[i]*σ[i])
    #         sum_of_weights += 1 / (σ[i]*σ[i])
    #     end
    #     wx = sum_of_values / sum_of_weights
    #
    #     @avx for i=1:n
    #         χ2 += (x[i] - wx) * (x[i] - wx) / (σ[i] * σ[i])
    #     end
    #     mswd = χ2 / (n-1)
    #     wσ = sqrt(1.0 / sum_of_weights)
    #     return wx, wσ, mswd
    # end
    export awmean

    """
    ```julia
    (wx, wσ, mswd) = gwmean(x, σ)
    ```
    Geochronologist's weighted mean, with "MSWD correction" to uncertainty,
    i.e., wσ is increased by a factor of sqrt(mswd)
    """
    function gwmean(x, σ)
        n = length(x)

        sum_of_values = sum_of_weights = χ2 = 0.0
        @inbounds @simd for i=1:n
            sum_of_values += x[i] / (σ[i]*σ[i])
            sum_of_weights += 1 / (σ[i]*σ[i])
        end
        wx = sum_of_values / sum_of_weights

        @inbounds @simd for i=1:n
            χ2 += (x[i] - wx) * (x[i] - wx) / (σ[i] * σ[i])
        end
        mswd = χ2 / (n-1)
        wσ = sqrt(mswd / sum_of_weights)
        return wx, wσ, mswd
    end
    # function gwmean(x::Array{<:Number}, σ::Array{<:Number})
    #     n = length(x)
    #     sum_of_values = sum_of_weights = χ2 = 0.0
    #     @avx for i=1:n
    #         sum_of_values += x[i] / (σ[i]*σ[i])
    #         sum_of_weights += 1 / (σ[i]*σ[i])
    #     end
    #     wx = sum_of_values / sum_of_weights
    #
    #     @avx for i=1:n
    #         χ2 += (x[i] - wx) * (x[i] - wx) / (σ[i] * σ[i])
    #     end
    #     mswd = χ2 / (n-1)
    #     wσ = sqrt(mswd / sum_of_weights)
    #     return wx, wσ, mswd
    # end
    export gwmean

    """
    ```julia
    MSWD(x, σ)
    ```
    Return the Mean Square of Weighted Deviates (AKA the reduced chi-squared
    statistic) of a dataset with values `x` and one-sigma uncertainties `σ`
    """
    function MSWD(x, σ)
        sum_of_values = sum_of_weights = χ2 = 0.0
        n = length(x)

        @inbounds @simd for i=1:n
            w = 1 / (σ[i]*σ[i])
            sum_of_values += w * x[i]
            sum_of_weights += w
        end
        wx = sum_of_values / sum_of_weights

        @inbounds @simd for i=1:n
            χ2 += (x[i] - wx) * (x[i] - wx) / (σ[i] * σ[i])
        end

        return χ2 / (n-1)
    end
    # function MSWD(x::Array{<:Number}, σ::Array{<:Number})
    #     sum_of_values = sum_of_weights = χ2 = 0.0
    #     n = length(x)
    #
    #     @avx for i=1:n
    #         w = 1 / (σ[i]*σ[i])
    #         sum_of_values += w * x[i]
    #         sum_of_weights += w
    #     end
    #     wx = sum_of_values / sum_of_weights
    #
    #     @avx for i=1:n
    #         χ2 += (x[i] - wx) * (x[i] - wx) / (σ[i] * σ[i])
    #     end
    #
    #     return χ2 / (n-1)
    # end
    export MSWD


## -- The York (1968) two-dimensional linear regression with x and y uncertainties
    # as commonly used in isochrons

    # Custom type to hold York fit resutls
    struct YorkFit{T}
        intercept::T
        intercept_sigma::T
        slope::T
        slope_sigma::T
        mswd::T
    end

    """
    ```julia
    yorkfit(x, σx, y, σy)
    ```
    Uses the York (1968) two-dimensional least-squares fit to calculate `a`, `b`,
    and uncertanties `σa`, `σb` for the equation `y = a + bx`, given `x`, `y` and
    uncertaintes `σx`, ``σy`.
    """
    function yorkfit(x, σx, y, σy; niterations=10)

        ## 1. Ordinary linear regression (to get a first estimate of slope and intercept)

        # Check for missing data
        t = (x.==x) .& (y.==y) .& (σx.==σx) .& (σy.==σy)
        x = x[t]
        y = y[t]
        σx = σx[t]
        σy = σy[t]

        # Calculate the ordinary least-squares fit
        # For the equation y=a+bx, m(1)=a, m(2)=b
        g = [ones(length(x)) x]
        m = (g'*g)\g'*y
        b = m[2]
        a = m[1]

        ## 2. Now, let's define parameters needed by the York fit

        # Weighting factors
        ωx = 1.0 ./ σx.^2
        ωy = 1.0 ./ σy.^2

        # terms that don't depend on a or b
        α = sqrt.(ωx .* ωy)

        x̄ = _vmean(x, :)
        ȳ = _vmean(y, :)
        r = sum((x .- x̄).*(y .- ȳ)) ./ (sqrt(sum((x .- x̄).^2)) * sqrt(sum((y .- ȳ).^2)))

        ## 3. Perform the York fit (must iterate)
        W = ωx.*ωy ./ (b^2*ωy + ωx - 2*b*r.*α)

        X̄ = sum(W.*x) / sum(W)
        Ȳ = sum(W.*y) / sum(W)

        U = x .- X̄
        V = y .- Ȳ

        sV = W.^2 .* V .* (U./ωy + b.*V./ωx - r.*V./α)
        sU = W.^2 .* U .* (U./ωy + b.*V./ωx - b.*r.*U./α)
        b = sum(sV) ./ sum(sU)

        a = Ȳ - b .* X̄
        for i = 2:niterations
            W .= ωx.*ωy ./ (b^2*ωy + ωx - 2*b*r.*α)

            X̄ = sum(W.*x) / sum(W)
            Ȳ = sum(W.*y) / sum(W)

            U .= x .- X̄
            V .= y .- Ȳ

            sV .= W.^2 .* V .* (U./ωy + b.*V./ωx - r.*V./α)
            sU .= W.^2 .* U .* (U./ωy + b.*V./ωx - b.*r.*U./α)
            b = sum(sV) ./ sum(sU)

            a = Ȳ - b .* X̄
        end

        ## 4. Calculate uncertainties and MSWD
        β = W .* (U./ωy + b.*V./ωx - (b.*U+V).*r./α)

        u = X̄ .+ β
        v = Ȳ .+ b.*β

        xm = sum(W.*u)./sum(W)
        ym = sum(W.*v)./sum(W)

        σb = sqrt(1.0 ./ sum(W .* (u .- xm).^2))
        σa = sqrt(1.0 ./ sum(W) + xm.^2 .* σb.^2)

        # MSWD (reduced chi-squared) of the fit
        mswd = 1.0 ./ length(x) .* sum( (y .- a.-b.* x).^2 ./ (σy.^2 + b.^2 .* σx.^2) )

        ## Results
        return YorkFit(a, σa, b, σb, mswd)
    end
    export yorkfit

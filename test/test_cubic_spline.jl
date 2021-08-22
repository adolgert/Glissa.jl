using Glissa
using Test
using SafeTestsets
using GLM
using DataFrames
using RCall


@safetestset "Cubic spline matches at points" begin
    x = collect(0:0.1:pi)
    f = cos.(x)
    cs = cubic_spline_flat_endpoints(x, f)
    for i1 in 1:length(x)
        @test abs(cs(x[i1]) - f[i1]) < 1e-9
    end
end


@safetestset "Cubic spline error reduces as dx reduces" begin
    Δ = (1e-5) * 2 .^ (0:14)
    maxerr = zeros(Float64, length(Δ))
    for j in 1:length(Δ)
        x = collect(0:(Δ[j]):(pi/2))
        f = cos.(x)
        cs = cubic_spline_flat_endpoints(x, f)

        xx = collect(0:(Δ[j]/10):(pi/2))
        mm = zero(Float64)
        for x1 in xx
            mm = max(mm, abs(cs(x1) - cos(x1)))
        end
        maxerr[j] = mm
    end
    df = DataFrame(logx = log.(Δ), logy = log.(maxerr))
    regress = lm(@formula(logy ~ logx), df)
    slope = coef(regress)[2]
    @test 0.7 < slope
    @test slope < 1.0
end


@safetestset "Cubic spline is continuous" begin
    Δ = 0.1
    x = collect(0:Δ:pi)
    f = sin.(x)
    cs = cubic_spline_flat_endpoints(x, f)
    for i1 in 1:(length(x) - 1)
        d = Δ / 10
        prev = 100
        # Assert values approach the expected value.
        for xl in range(x[i1] + Δ / 10, length = 5, stop = x[i1] + 1e-9)
            v = abs(sin(x[i1]) - cs(xl))
            @test v <= prev
            prev = v
        end
        # Assert they get close-enough for the final value.
        @test prev < 1e-5
    end
    for i1 in 2:length(x)
        d = Δ / 10
        prev = 100
        # Look at approach from below.
        for xl in range(x[i1] - Δ / 10, length = 5, stop = x[i1] - 1e-9)
            v = abs(sin(x[i1]) - cs(xl))
            @test v <= prev
            prev = v
        end
        @test prev < 1e-5
    end
end



@safetestset "Cubic spline is continuous for uneven x" begin
    Δ = 0.1
    x = vcat([0.0, 0.1, 0.2, 0.4, 0.5, 1], 2:5)
    f = sin.(x)
    cs = cubic_spline_flat_endpoints(x, f)
    for i1 in 1:(length(x) - 1)
        d = Δ / 10
        prev = 100
        # Assert values approach the expected value.
        for xl in range(x[i1] + Δ / 10, length = 5, stop = x[i1] + 1e-9)
            v = abs(sin(x[i1]) - cs(xl))
            @test v <= prev
            prev = v
        end
        # Assert they get close-enough for the final value.
        @test prev < 1e-5
    end
    for i1 in 2:length(x)
        d = Δ / 10
        prev = 100
        # Look at approach from below.
        for xl in range(x[i1] - Δ / 10, length = 5, stop = x[i1] - 1e-9)
            v = abs(sin(x[i1]) - cs(xl))
            @test v <= prev
            prev = v
        end
        @test prev < 1e-5
    end
end


@safetestset "Cubic spline overshoot" begin
    # RPN 15A data from Fritsch-Carlson
    x = [7.99, 8.09, 8.19, 8.7, 9.2, 10, 12, 15, 20]
    f = [0, 2.76429e-5, 4.37498e-2, 0.169183, 0.469428, 0.943740, 0.998636, 0.999919, 0.999994]
    cs = cubic_spline(x, f, [0.0, 0.0])
    xx = 8.08:0.0001:8.1
    plot(xx, cs.(xx))
    ϵ = 1e-8
    for i in 2:(length(x) - 1)
        left = (cs(x[i] - ϵ) - cs(x[i] - 2ϵ)) / ϵ
        right = (cs(x[i] + 2ϵ) - cs(x[i] + ϵ))/ϵ
        d = abs(left - right)
        @test d < 1e-6
    end
    ϵ = 1e-5
    second_derivative = zeros(Float64, length(x) - 2)
    for i in 2:(length(x) - 1)
        left = (cs(x[i] - 3ϵ) - 2*cs(x[i] - 2ϵ) + cs(x[i] - ϵ)) / ϵ^2
        right = (cs(x[i] + ϵ) - 2*cs(x[i] + 2ϵ) + cs(x[i] + 3ϵ)) / ϵ^2
        second_derivative[i - 1] = abs(left - right)
    end
    @test second_derivative[1] > 5 * second_derivative[5]
end


@safetestset "Hyman criterion matches splinefun in R" begin
    # this hyman filter comes from the R splinefun source code.
    # x is the knots. y is values. b is gradient.
    R"""hyman_filter <- function(x, y, b)
    {
        n <- length(x)
        ss <- diff(y) / diff(x)
        S0 <- c(ss[1L], ss)
        S1 <- c(ss, ss[n-1L])
        t1 <- pmin(abs(S0), abs(S1))
        sig <- b
        ind <- S0*S1 > 0
        sig[ind] <- S1[ind]
        ind <- sig >= 0
        if(sum(ind)) b[ind] <- pmin(pmax(0, b[ind]), 3*t1[ind])
        ind <- !ind
        if(sum(ind)) b[ind] <- pmax(pmin(0, b[ind]), -3*t1[ind])
        b
    }
    """

    x = 0:0.1:π
    y = sin.(x)
    yp = cos.(x)
    hf = collect(R"""hyman_filter($x, $y, $yp)""")
    Glissa.project_to_monotonicity!(x, y, yp)
    yp .- hf
    @test maximum(abs.(yp .- hf)) < 1e-12
end

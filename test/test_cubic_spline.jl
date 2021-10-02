using DataFrames
using GLM
using RCall
using Test


@testset "Cubic spline matches at points" begin
    using Glissa
    x = collect(0:0.1:pi)
    f = cos.(x)
    cs = cubic_spline(x, f, [0, 0]; slope=FreeSlope())
    for i1 in 1:length(x)
        @test abs(cs(x[i1]) - f[i1]) < 1e-9
    end
end

function points_slope(logx, logy)
    df = DataFrame(logx = logx, logy = logy)
    model = @formula(logy ~ logx)
    regress = lm(model, df)
    slope = coef(regress)[2]
    slope
end

@testset "Cubic spline error reduces as dx reduces" begin
    using Glissa
    Δ = (1e-5) * 2 .^ (0:14)
    maxerr = zeros(Float64, length(Δ))
    for j in 1:length(Δ)
        x = collect(0:(Δ[j]):(pi/2))
        f = cos.(x)
        cs = cubic_spline(x, f, [0, 0]; slope=FreeSlope())

        xx = collect(0:(Δ[j]/10):(pi/2))
        mm = zero(Float64)
        for x1 in xx
            mm = max(mm, abs(cs(x1) - cos(x1)))
        end
        maxerr[j] = mm
    end
    slope = points_slope(log.(Δ), log.(maxerr))
    @test 0.7 < slope
    @test slope < 1.0
end


@testset "Cubic spline is continuous" begin
    using Glissa
    Δ = 0.1
    x = collect(0:Δ:pi)
    f = sin.(x)
    cs = cubic_spline(x, f, [0, 0]; slope=FreeSlope())
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


@testset "Cubic spline is continuous for uneven x" begin
    using Glissa
    Δ = 0.1
    x = vcat([0.0, 0.1, 0.2, 0.4, 0.5, 1], 2:5)
    f = sin.(x)
    cs = cubic_spline(x, f, [0, 0]; slope=FreeSlope())
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


@testset "Cubic spline overshoot" begin
    using Glissa
    # RPN 15A data from Fritsch-Carlson
    x = [7.99, 8.09, 8.19, 8.7, 9.2, 10, 12, 15, 20]
    f = [0, 2.76429e-5, 4.37498e-2, 0.169183, 0.469428, 0.943740, 0.998636, 0.999919, 0.999994]
    cs = cubic_spline(x, f, [0.0, 0.0]; slope = Monotonic())
    xx = 8.08:0.0001:8.1
    # plot(xx, cs.(xx))
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


@testset "Hyman criterion matches splinefun in R" begin
    using Glissa
    using RCall
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
    x = collect(0:0.1:π)
    y = sin.(x)
    yp = cos.(x)
    hf = collect(R"""hyman_filter($x, $y, $yp)""")
    Glissa.project_to_monotonicity!(x, y, yp)
    @test maximum(abs.(yp .- hf)) < 1e-12

    # This example is constructed by hand to need monotonic smoothing.
    x = [7.99, 8.09, 8.19, 8.7, 9.2, 10, 12, 15, 20]
    f = [0, 2.76429e-5, 4.37498e-2, 0.169183, 0.469428, 0.943740, 0.998636, 0.999919, 0.999994]
    fp = similar(f)
    Glissa.global_derivatives!(x, f, fp)
    hf = collect(R"""hyman_filter($x, $f, $fp)""")
    # Assert that the uncorrected slopes are different from corrected ones.
    @test maximum(abs.(hf - fp)) > 1e-10
    Glissa.project_to_monotonicity!(x, f, fp)
    # After correction, they agree with the R version.
    @test maximum(abs.(fp .- hf)) < 1e-12

    R"""spl_coef_conv <- function(x, y, b)
    {
        n <- length(x)
        h <- diff(x); y <- -diff(y)
        b0 <- b[-n]; b1 <- b[-1L]
        cc <- -(3*y + (2*b0 + b1)*h) / h^2
        c1 <- (3*y[n-1L] + (b0[n-1L] + 2*b1[n-1L])*h[n-1L]) / h[n-1L]^2
        c <- c(cc, c1)
        dd <- (2*y/h + b0 + b1) / h^2
        d <- c(dd, dd[n-1L])
        list(c, d)
    }
    """
    cd = R"""spl_coef_conv($x, $f, $fp)"""
    c = collect(cd[1])
    d = collect(cd[2])
    C = zeros(Float64, 4, length(x) - 1)
    Glissa.hyman_coefficients!(x, f, fp, C)
    # The R functions define vectors that are 1 longer.
    @test maximum(abs.(C[3, :] - c[1:(end-1)])) < 1e-12
    @test maximum(abs.(C[4, :] - d[1:(end-1)])) < 1e-12
end


@testset "Cubic spline is compatible with forward diff" begin
    using Glissa
    using ForwardDiff
    x = [7.99, 8.09, 8.19, 8.7, 9.2, 10, 12, 15, 20]
    x2 = collect(8:0.01:10)
    y = [0, 2.76429e-5, 4.37498e-2, 0.169183, 0.469428, 0.943740, 0.998636, 0.999919, 0.999994]

    function g(f::AbstractVector{T}) where {T}
        cs = Glissa.cubic_spline(x, f, zeros(T, 2); slope = Monotonic())
        cs.(x2) 
    end

    interp = g(y)
    jacob = ForwardDiff.jacobian(g, y)
    @test size(jacob) == (201, length(x))
    @test jacob[1, 1] > 0
    @test abs(jacob[1, 9]) < 1e-7
end


@testset "Cubic spline is compatible with reverse diff" begin
    using Glissa
    using ReverseDiff
    x = [7.99, 8.09, 8.19, 8.7, 9.2, 10, 12, 15, 20]
    x2 = collect(8:0.01:10)
    f = [0, 2.76429e-5, 4.37498e-2, 0.169183, 0.469428, 0.943740, 0.998636, 0.999919, 0.999994]
    function g(f::AbstractVector{T}, x1, x2) where {T}
        cs = Glissa.cubic_spline(x1, f, zeros(T, 2); slope = Monotonic())
        cs.(x2) 
    end
    interp = g(f, x, x2)
    f_tape = ReverseDiff.GradientTape(g, (f, x, x2))
    compiled_f_tape = ReverseDiff.compile(f_tape)
    inputs = (f, x, x2)
    cfg = ReverseDiff.GradientConfig(inputs)
    results = similar(interp)
    # This doesn't work, and I'd like it to work.
    # ReverseDiff.gradient!(results, compiled_f_tape, inputs)
end


@testset "Cubic spline integrates a flat line" begin
    using Glissa
    using Random
    using Distributions
    x = Float64[0, 1, 3, 5]
    y = Float64[1, 1, 1, 1]
    cs = Glissa.cubic_spline(x, y, [0.0, 0.0]; slope = Monotonic())
    @test abs(cs(2.5) - 1) < 1e-6
    exact(z) = z
    rng = MersenneTwister(9234724)
    for j in 1:100
        a = rand(rng, Uniform(0, 5))
        b = rand(rng, Uniform(a, 5))
        found = Glissa.integrate(cs, a, b)
        @test abs(found - exact(b) + exact(a)) < 1e-9
    end
end

@testset "Cubic spline integrates a sloped line" begin
    using Glissa
    using Random
    using Distributions
    x = Float64[0, 1, 3, 5]
    y = Float64[0, 1, 3, 5]
    cs = Glissa.cubic_spline(x, y, [1.0, 1.0]; slope = Monotonic())
    @test abs(cs(2.5) - 2.5) < 1e-6
    exact(z) = 0.5 * z^2
    rng = MersenneTwister(9292734)
    for j in 1:100
        a = rand(rng, Uniform(0, 5))
        b = rand(rng, Uniform(a, 5))
        found = Glissa.integrate(cs, a, b)
        @test abs(found - exact(b) + exact(a)) < 1e-9
    end
end


@testset "Cubic spline integrates a quadratic line" begin
    using Glissa
    using Random
    using Distributions
    Δ = 0.1
    x = collect(0:Δ:5)
    y = x.^2
    cs = Glissa.cubic_spline(x, y, [0.0, 10.0]; slope = Monotonic())
    for q in [0.3, 1.2, 4.7, 4.99]
        @test abs(cs(q) - q^2) < Δ^3
    end
    exact(z) = z^3 / 3
    rng = MersenneTwister(89192734)
    for j in 1:100
        a = rand(rng, Uniform(0, 5))
        b = rand(rng, Uniform(a, 5))
        found = Glissa.integrate(cs, a, b)
        @test abs(found - exact(b) + exact(a)) < 8*Δ^3
    end
end


@testset "Cubic spline integrates a cosine" begin
    using Glissa
    using Random
    using Distributions
    Δ = π / 30
    x = collect(0:Δ:π)
    y = cos.(x)
    cs = Glissa.cubic_spline(x, y, [0.0, 0.0]; slope = Monotonic())
    for q in [0.3, 1.2, 1.7, 2.99]
        @test abs(cs(q) - cos(q)) < Δ^3
    end
    exact(z) = sin(z)
    rng = MersenneTwister(89192734)
    for j in 1:100
        a = rand(rng, Uniform(0, 3))
        b = rand(rng, Uniform(a, 3))
        found = Glissa.integrate(cs, a, b)
        @test abs(found - exact(b) + exact(a)) < Δ^3
    end
end

@testset "Cubic spline derivative" begin
    using Glissa
    cs = CubicSpline{Int,Float64}([0, 1, 2],
        Float64[
            1 2;
            3 7;
            4 0;
            1 -1
        ]
    )
    csprime = Glissa.derivative(cs)
    @test size(csprime.c) == (3, 2)
    @test csprime.c ≈ Float64[
        3 7;
        8 0;
        3 -3
    ]
end

@testset "Cubic spline multiplication" begin
    using Glissa
    cs1 = CubicSpline{Int,Float64}([0, 1, 2],
        Float64[
            1 2;
            3 7;
            4 0
        ]
    )
    using Glissa
    cs2 = CubicSpline{Int,Float64}([0, 1, 2],
        Float64[
            2 2;
            1 -2
        ]
    )
    cs3 = cs1 * cs2
    @test cs3.c ≈ Float64[
        2 4;
        7 10;
        11 -14;
        4 0
    ]
end
using Random
using Test
using Revise
using Glissa

# There are three implementations of divided differences.
# 1. By a ratio of determinants: divdiff
# 2. By recursion: divided_difference
# 3. Explicit: divided_difference_explicit. (distinct axis knots only)
#
# We want to see that these all agree.

@testset "unique_knots Float64 all methods agree" begin
    rng = Random.MersenneTwister(98742342)
    n = 10
    for i in 1:1000
        axis = sort(10 * rand(rng, n))
        ukf(x) = cos(x) * (1 - 0.1 * x^2)
        # Pick a random sub-range of the axis on which to compute the divided difference.
        sub = rand(rng, 1:4)
        len = rand(rng, 0:5)
        subr = sub:(sub + len)
        dd_det = Glissa.divdiff(view(axis, subr), ukf)
        dd = Glissa.DividedDifference{Float64}()
        dd_rec = Glissa.divided_difference(subr.start, subr.stop, axis, ukf, dd)
        dd_exp = Glissa.divided_difference_explicit(subr.start, subr.stop, axis, ukf)
        # Check the three versions are equal.
        # Why is this eps so large? It's a relative error, and the calculation using
        # determinants is prone to disagreement. It's not very stable.
        eps = 1e-6
        abs1 = abs((dd_det - dd_rec) / dd_det)
        abs2 = abs((dd_det - dd_exp) / dd_det)
        if abs1 ≥ eps || abs2 ≥ eps
            @show i, subr, dd_det, dd_rec, dd_exp, abs1, abs2
            @test abs1 < eps
            @test abs2 < eps
        end
    end
end

@testset "unique_knots bigfloat comparison" begin
    rng = Random.MersenneTwister(728432342)
    float_precision = precision(BigFloat)
    setprecision(BigFloat, 256)
    maxdiff = zeros(BigFloat, 3)
    n = 10
    for i in 1:1000
        axis = sort(10 * rand(rng, Float64, n))
        axis_big = convert(Vector{BigFloat}, axis)
        ukf(x) = sin(x) * (3//2 - 1//10 * x^3)
        # Pick a random sub-range of the axis on which to compute the divided difference.
        sub = rand(rng, 1:4)
        len = rand(rng, 0:5)
        subr = sub:(sub + len)
        dd_det = Glissa.divdiff(view(axis, subr), ukf)
        ddmemo = Glissa.DividedDifference{Float64}()
        dd_rec = Glissa.divided_difference(subr.start, subr.stop, axis, ukf, ddmemo)
        dd_exp = Glissa.divided_difference_explicit(subr.start, subr.stop, axis, ukf)
        dd_big = Glissa.divided_difference_explicit(subr.start, subr.stop, axis_big, ukf)
        eps = 1e-7
        # I checked all three with BigFloat, and they agree to 1e-59, or something,
        # so the dd_big is the true number.
        absdd = [
            abs((dd_det - dd_big) / dd_big)
            abs((dd_rec - dd_big) / dd_big)
            abs((dd_exp - dd_big) / dd_big)
        ]
        for j in 1:3
            if absdd[j] > maxdiff[j]
                maxdiff[j] = absdd[j]
            end
        end
        # if abs1 ≥ eps || abs2 ≥ eps
        #     @test abs1 < eps
        #     @test abs2 < eps
        # end
    end
    @show maxdiff
    setprecision(BigFloat, float_precision)
    @test maxdiff[1] < 1e-8
    @test maxdiff[3] < maxdiff[1]  # ~1e-10
    @test maxdiff[2] < maxdiff[3]  # ~1e-10
end

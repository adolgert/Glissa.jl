using Test
using Glissa
using Random
using Distributions
using LinearAlgebra


@testset "recursive Q B-spline has correct positivity" begin
    rng = Random.MersenneTwister(23429235)
    for trial_idx in 1:1000
        axis, uniques, multiples, order = Glissa.generate_random_polyspline_specification(rng, Float64)
        bspline_cnt = sum(multiples) - order
        index = rand(rng, 1:bspline_cnt)
        # Ensure it's 0 to the left of the B-spline
        if index > 1
            xleft = axis[1]
            bleft = Glissa.bsplineq_recursive(axis, order, index, xleft)
            @test bleft == 0
        end
        # Ensure it's 0 to the right of the B-spline
        if index + order < length(axis)
            xright = axis[end]
            bright = Glissa.bsplineq_recursive(axis, order, index, xright)
            @test bright == 0
        end
        # Ensure it's positive within the B-spline
        # The points are internal to the B-spline. Excludes endpoints.
        pt_cnt = 10
        dy = (axis[index + order] - axis[index]) / (pt_cnt + 1)
        test_points = [axis[index] + i * dy for i in 1:pt_cnt]
        for pt in test_points
            pos1 = Glissa.bsplineq_recursive(axis, order, index, pt)
            @test pos1 > 0
        end
    end
end

@testset "Qspline has correct first moment" begin
    rng = Random.MersenneTwister(23429235)
    for trial_idx in 1:10
        axis, uniques, multiples, order = Glissa.generate_random_polyspline_specification(rng, Float64)
        bspline_cnt = sum(multiples) - order
        bspline_index = rand(rng, 1:bspline_cnt)
        a = axis[bspline_index]
        b = axis[bspline_index + order]
        dx = (b - a) / 100000
        total = 0.0
        x = a
        while x ≤ b
            pos1 = Glissa.bsplineq_recursive(axis, order, bspline_index, x)
            total += dx * pos1
            x += dx
        end
        @test abs(total - 1 / order) * order < 1e-4
    end
end


# Create test points inside every interval of a nonuniform axis.
# The axis has no repeated values.
function generate_test_points(axis, pt_cnt)
    pts = zeros(eltype(axis), (length(axis) - 1) * pt_cnt)
    pt_idx = 1
    for interval_idx in 1:(length(axis) - 1)
        dy = (axis[interval_idx + 1] - axis[interval_idx]) / (pt_cnt + 1)
        for j in 1:pt_cnt
            pts[pt_idx] = axis[interval_idx] + dy * j
            pt_idx += 1
        end
    end
    pts
end


# See the documentation for "testspaces.md". It explains this.
@testset "recursive Q B-spline is a basis for a polyspline space" begin
    rng = Random.MersenneTwister(23429235)
    for trial_idx = 1:100
        T = Float64
        axis, uniques, multiples, order = Glissa.generate_random_polyspline_specification(rng, T)
        uniques_cnt = length(uniques)
        poly_coeff_cnt = (uniques_cnt - 1) * order
        A = Glissa.polyspline_constraints!(uniques, multiples, order)
        # b is the right-hand-side of the Ac=b equation, and it's all zeroes for these constraints.
        b = zeros(T, size(A, 1))
        @assert size(A, 2) == poly_coeff_cnt
        Astar = LinearAlgebra.pinv(A)

        # Create a random guess for the coefficients of a polynomial spline.
        w = 2*rand(rng, T, poly_coeff_cnt) .- 1
        # Project out the disallowed parts of those coefficients to find a random polyspline c.
        c = Astar * b + (Diagonal(ones(poly_coeff_cnt)) - Astar * A) * w
        polyspline = Glissa.PolySpline{T,T}(uniques, reshape(c, (order, uniques_cnt - 1)))

        pts = generate_test_points(uniques, order + 1)
        ytrue = similar(pts)
        for genidx in 1:length(pts)
            ytrue[genidx] = polyspline(pts[genidx])
        end

        bspline_cnt = sum(multiples) - order
        M = zeros(T, length(pts), bspline_cnt)
        for bidx in 1:bspline_cnt
            for trial_idx in 1:length(pts)
                M[trial_idx, bidx] = Glissa.bsplineq_recursive(axis, order, bidx, pts[trial_idx])
            end
        end

        bspline_coeffs = LinearAlgebra.pinv(M) * ytrue
        err = zeros(T, length(ytrue))
        for check_idx in 1:length(pts)
            yt = zero(T)
            for bidx in 1:bspline_cnt
                yt += bspline_coeffs[bidx] * Glissa.bsplineq_recursive(axis, order, bidx, pts[check_idx])
            end
            err[check_idx] = yt - ytrue[check_idx]
        end
        # @show err
        @test maximum(abs.(err)) < 1e-12
    end
end

# There are lots of B-spline implementations, so let's see if they all agree.
using Random
using Distributions
using Test
using Revise
using Glissa


@testset "Value agrees when comparing implementations" begin
    # Compares recursive b-splines with schumaker's algorithm 5.5.
    # Algorithm 5.5 expects padding of identical points near edges. That padding means
    # it produces extra B-splines near edges.
    rng = Random.MersenneTwister(3224923)
    T = Float64
    largest_err = zero(T)
    for trial_idx in 1:100
        # [1, 1, 3, 1] order 5. Causing trouble.
        axis, uniques, multiples, order = Glissa.generate_random_polyspline_specification(rng, T)
        if order == 1  # Algorithm 5.5 can't handle order=1.
            continue
        end
        # Schumaker's algorithms assume there is padding on the sides of repeated knots.
        multiples[1] = order
        multiples[end] = order
        axis = Glissa.axis_multiples_to_repeats(uniques, multiples)
        bspline_cnt = Glissa.bspline_count(axis, order)
        interval_cnt = length(uniques) - 1
        # This is an array of which bsplines are in which interval.
        cover = Glissa.bspline_indices_in_interval(multiples, order)
        ans_rec = zeros(T, order)
        ans_iter = zeros(T, order + 1) # The +1 because this needs extra buffer space.
        for interval_idx in 1:interval_cnt
            a, b = (uniques[interval_idx], uniques[interval_idx + 1])
            n = order + 1  # Need some points to distinguish polynomials.
            dx = (b - a) / (n + 1)
            for ix in 1:n
                x = a + dx * ix
                ans_rec = zeros(T, order)

                # Recursive
                ans_rec .= 0
                for ci in 1:order
                    if cover[ci, interval_idx] > 0
                        bidx = cover[ci, interval_idx]
                        ans_rec[ci] = Glissa.bsplinen_recursive(axis, order, bidx, x)
                    end
                end
                l = searchsortedlast(axis, x)
                # Use the padded axis here.
                Glissa.generate_normalized_bsplines!(ans_iter, axis, l, order, x)
                @show ans_rec, ans_iter
                err = maximum(abs.((ans_rec - ans_iter[1:order]) ./ ans_rec))
                largest_err = max(err, largest_err)
                @test err < 1e-13
            end
        end
    end
    @show largest_err
end

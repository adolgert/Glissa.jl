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
    for trial_idx in 1:100
        # [1, 1, 3, 1] order 5. Causing trouble.
        axis, uniques, multiples, order = Glissa.generate_random_polyspline_specification(rng, T)
        if order == 1  # Algorithm 5.5 can't handle order=1.
            continue
        end
        bspline_cnt = Glissa.bspline_count(axis, order)
        interval_cnt = length(uniques) - 1
        # This is an array of which bsplines are in which interval.
        cover = Glissa.bspline_indices_in_interval(multiples, order)
        ans_rec = zeros(T, order)
        ans_iter = zeros(T, order + 1) # The +1 because this needs extra buffer space.
        augmented_mult = copy(multiples)
        # Schumaker's algorithms assume there is padding on the sides of repeated knots.
        augmented_mult[1] = order
        augmented_mult[end] = order
        augmented = Glissa.axis_multiples_to_repeats(uniques, augmented_mult)
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
                l = searchsortedlast(augmented, x)
                # Use the padded axis here.
                Glissa.generate_normalized_bsplines!(ans_iter, augmented, l, order, x)
                @show ans_rec, ans_iter
                for eg_idx in 1:order
                    if cover[eg_idx, interval_idx] > 0
                        ans_rec1 = ans_rec[eg_idx]
                        rec1_found = false
                        for check_iter in 1:order
                            if abs((ans_iter[check_iter] - ans_rec1) / ans_rec1) < 1e-6
                                rec1_found = true
                            end
                        end
                        if !rec1_found
                            @show order, multiples, uniques
                            @show cover
                            @show x, ans_rec, ans_iter
                            @test rec1_found
                        end
                    end
                end
            end
        end
    end
end

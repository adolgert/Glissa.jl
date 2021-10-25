using Test
using Glissa
using Random
using Distributions


@testset "recursive Q B-spline has correct positivity" begin
    rng = Random.MersenneTwister(23429235)
    for i in 1:1000
        order = rand(rng, 1:5)
        index = rand(rng, 1:5)
        axis_cnt = rand(rng, (index + order):(index + order + 3))
        uniques = sort(rand(rng, Float64, axis_cnt))
        mult_type = rand(rng, [:Single, :Same, :Random])
        multiples = ones(Int, length(uniques))
        if mult_type == :Same && order > 1
            multiples .= rand(rng, 2:order)
        elseif mult_type == :Random
            p = vcat(ones(Int, order - 1), [0.1]) # discount the disconnected polynomial
            multiples .= rand(rng, Categorical(p/sum(p)), length(multiples))
        end
        axis = Glissa.axis_multiples_to_repeats(uniques, multiples)
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


@testset "recursive Q B-spline is a space" begin
    rng = Random.MersenneTwister(23429235)
    for i in 1:1000
        order = rand(rng, 1:5)
        index = rand(rng, 1:5)
        axis_cnt = rand(rng, (index + order):(index + order + 3))
        uniques = sort(rand(rng, Float64, axis_cnt))
        mult_type = rand(rng, [:Single, :Same, :Random])
        multiples = ones(Int, length(uniques))
        if mult_type == :Same && order > 1
            multiples .= rand(rng, 2:order)
        elseif mult_type == :Random
            p = vcat(ones(Int, order - 1), [0.1]) # discount the disconnected polynomial
            multiples .= rand(rng, Categorical(p/sum(p)), length(multiples))
        end
        axis = Glissa.axis_multiples_to_repeats(uniques, multiples)
        bspline_cnt = sum(multiples) - order
        for bidx in 1:bspline_cnt
            for axidx in 1:(length(axis) - 1)
                q = Glissa.bsplineq_recursive(axis, order, index, axis[bidx])
            end
        end
    end
end

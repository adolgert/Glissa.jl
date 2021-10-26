using Test
using Glissa

@testset "Explicit examples of switching from repeats to multiples and back" begin
    examples = [
        ([1.0, 1.0, 2.0, 2.0, 3.0, 3.0], [1.0, 2.0, 3.0], [2, 2, 2]),
        (Float16[1.0, 2.0, 3.0, 4.0], Float16[1.0, 2.0, 3.0, 4.0], [1, 1, 1, 1]),
        (Float64[], Float64[], Int[]),
        ([1.1], [1.1], [1]),
        (BigFloat[1.2, 1.2], BigFloat[1.2], [2]),
        ([1.1, 1.2, 1.2, 1.2, 1.2], [1.1, 1.2], [1, 4])
    ]
    for (repeated, unique, multiples) in examples
        un1, mult1 = Glissa.axis_repeats_to_multiples(repeated)
        @test un1 == unique
        @test mult1 == multiples
        @test eltype(un1) == eltype(repeated)

        rep1 = Glissa.axis_multiples_to_repeats(unique, multiples)
        @test rep1 == repeated
        @test eltype(rep1) == eltype(unique)
    end
end


@testset "bspline_indices_in_interval" begin
    mu = [2, 1, 3, 2, 2, 1]
    order = 3
    mi = Glissa.MultIndex(mu)
    mi[4]
    
    cover = bspline_indices_in_interval(mu, order)
end

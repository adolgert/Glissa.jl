using Glissa
using Test
using Random
using Distributions

@testset "the right indices are found" begin
    rng = Random.MersenneTwister(90832742)
    for trial_idx in 1:10
        arr_cnt = 10
        arr = sort(rand(rng, Uniform(-5, 5), arr_cnt))
        for subtrial_idx in 1:100
            x = rand(rng, Uniform(-5, 5))
            l = rand(rng, 1:(arr_cnt - 1))
            a = searchsortedlast(arr, x)
            b = Glissa.find_left_index(arr, x, l)
            @test  a == b
        end
    end
end

# Schumaker's math uses multiplicity to describe coincident knots,
# so m=1 for a single knot, and m=2 for a double-knot, up to m=order.
# This converts an axis with multiplicity to an axis that repeats identical values
# where knots are multiple, because Schumaker's algorithms expect this format.
function multiplicity_to_repetition(axis::AbstractVector{T}, m) where {T}
    y = zeros(T, sum(m))
    vertex_idx = 1
    multiplicity_idx = 1
    y_idx = 1
    while vertex_idx ≤ length(axis)
        y[y_idx] = axis[vertex_idx]
        y_idx += 1
        multiplicity_idx += 1
        if multiplicity_idx > m[vertex_idx]
            vertex_idx += 1
            multiplicity_idx = 1
        end
    end
    println("$(vertex_idx) $(length(axis)) $y")
    @assert vertex_idx == length(axis) + 1
    y
end

@testset "Algorithm 5.5 matches our explicit version" begin

    rng = Random.MersenneTwister(79234298)
    T = Float64
    order = 3 # order
    x = 5
    axis = convert(Vector{T}, 0:1:10) # axis
    multiplicity = ones(Int, length(axis))
    multiplicity[1] = order
    multiplicity[end] = order
    y = multiplicity_to_repetition(axis, multiplicity)

    # Use the explicit solution method.
    splines = Vector{PolyBSpline{T,T}}(undef, sum(multiplicity) - order)
    for i in 1:(sum(multiplicity) - order)
        ((l, r), m) = Glissa.reduce_axis(multiplicity, order, i)
        coeffs = zeros(T, order, r - l)
        println("$m l:r $(l):$(r)")
        Glissa.bspline_by_matrix!(axis[l:r], coeffs, m, order, :N)
        splines[i] = PolyBSpline{T,T}(axis[l:r], coeffs, UnitRange{Int}(l, r))
    end

    # Use Schumaker's method
    maxdiff = 0.0
    for trial_idx in 1:20
        x = rand(rng, Uniform(axis[1], axis[end]))
        l = searchsortedlast(y, x)
        N = zeros(T, order + 1)
        Glissa.generate_normalized_bsplines!(N, y, l, order, x)
        i = searchsortedlast(axis, x)
        println(N)
        match_idx = 1
        for j in 1:length(splines)
            if i ≥ splines[j].bounds.start && i < splines[j].bounds.stop
                value = splines[j](x)
                println("$j $(splines[j].bounds) $value $(match_idx)")
                maxdiff = max(maxdiff, abs(N[match_idx] - value))
                match_idx += 1
            end
        end
    end
    @test maxdiff < 1e-15

end

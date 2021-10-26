function axis_repeats_to_multiples(axis)
    uniques = copy(axis)
    multiplicity = zeros(Int, length(axis))
    ucnt = 0
    last = zero(eltype(axis))
    for uidx in 1:length(axis)
        if ucnt == 0 || axis[uidx] != last
            ucnt += 1
            uniques[ucnt] = axis[uidx]
            multiplicity[ucnt] = 1
            last = axis[uidx]
        else
            multiplicity[ucnt] += 1
        end
    end
    @assert sum(multiplicity) == length(axis)
    (uniques[1:ucnt], multiplicity[1:ucnt])
end


function axis_multiples_to_repeats(uniques, multiplicity)
    axis = zeros(eltype(uniques), sum(multiplicity))
    ax_idx = 1
    for uidx in 1:length(uniques)
        axis[ax_idx:(ax_idx + multiplicity[uidx] - 1)] .= uniques[uidx]
        ax_idx += multiplicity[uidx]
    end
    @assert ax_idx == length(axis) + 1
    axis
end


# Create a polyspline axis.
# This polyspline will have sum(multiples) - order Bsplines.
function generate_random_polyspline_specification(rng, T)
    order = rand(rng, 1:5)
    multiples_cnt = rand(rng, (order + 1):(2*order + 1))
    mult_type = rand(rng, [:Single, :Same, :Random])
    if mult_type == :Same && order > 1
        println("same")
        m_all = rand(rng, 2:order)
        m_cnt = div(multiples_cnt - 1, m_all) + 1
        multiples = ones(Int, m_cnt) * m_all
        multiples_cnt = length(multiples)
        @show mult_type, multiples_cnt, multiples
    elseif mult_type == :Random && order > 1
        println("random")
        # discount the disconnected polynomial
        p = vcat(ones(Int, order - 1), [0.1])
        dist = Categorical(p / sum(p))
        multiples = Int[]
        remain = multiples_cnt
        while remain > 0
            madd = rand(rng, dist)
            if madd > remain
                madd = remain
            end
            push!(multiples, madd)
            remain -= madd
        end
        @show mult_type, multiples_cnt, multiples
    else
        println("single")
        multiples = ones(Int, multiples_cnt)
        @show mult_type, multiples_cnt, multiples
    end
    println("uniques")
    uniques = sort(rand(rng, T, length(multiples)))
    axis = Glissa.axis_multiples_to_repeats(uniques, multiples)
    @assert length(axis) == sum(multiples)
    @assert length(uniques) == length(multiples)
    (axis, uniques, multiples, order)
end

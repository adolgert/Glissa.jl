using Random
using Distributions

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


bspline_count(axis::AbstractVector{T}, order) where {T <: Real} = length(axis) - order
bspline_count(multiples::AbstractVector{T}, order) where {T <: Integer} = sum(multiples) - order


struct MultIndex
    multiples::Vector{Int}
end

function Base.getindex(mi::MultIndex, i::Int)
    cm = cumsum(mi.multiples)
    j = searchsortedfirst(cm, i)
    if j ≤ length(mi.multiples)
        mult_idx = mi.multiples[j] - cm[j] + i
        (j, mult_idx)
    else
        (length(mi.multiples) + 1, 1)
    end
end


function getget(mi::MultIndex, i::Int)
    interval_cnt = length(mi.multiples) - 1
    interval_idx = 1
    multiple_idx = 1
    bidx = 1
    while interval_idx ≤ interval_cnt
        if bidx == i
            return (interval_idx, multiple_idx)
        end
        bidx += 1
        multiple_idx += 1
        if multiple_idx > mi.multiples[interval_idx]
            multiple_idx = 1
            interval_idx += 1
        end
    end
    return (interval_cnt + 1, 1)
end

function less_than_cover(i, j)
    if i == 0 && j == 0
        false
    elseif i == 0
        false
    elseif j == 0
        true
    else
        i < j
    end
end


function bspline_indices_in_interval(multiples, order)
    # Each spline covers between τ[i] and τ[i+order]
    mi = MultIndex(multiples)
    interval_cnt = length(multiples) - 1
    bspline_cnt = bspline_count(multiples, order)
    cover = zeros(Int, (order, interval_cnt))
    for bidx in 1:bspline_cnt
        start_axis_idx, start_mult_index = mi[bidx]
        finish_axis_idx, finish_mult_index = mi[bidx + order]
        start_interval_idx = start_axis_idx
        finish_interval_idx = finish_axis_idx - 1
        if finish_interval_idx > interval_cnt
            error("too few intervals $bidx $(bspline_cnt) $(finish_interval_idx)")
        end
        for iidx in start_interval_idx:finish_interval_idx
            slot = order
            while slot > 0 && cover[slot, iidx] != 0
                slot -= 1
            end
            if slot == 0
                error("more splines than the order, $bidx $iidx $cover")
            end
            cover[slot, iidx] = bidx
        end
    end
    # This is to match the order that Schumaker uses.
    sort(cover, dims=1, rev = true, lt=less_than_cover)
end


# Create a polyspline axis.
# This polyspline will have sum(multiples) - order Bsplines.
function generate_random_polyspline_specification(rng, T)
    order = rand(rng, 1:5)
    multiples_cnt = rand(rng, (order + 1):(2*order + 1))
    mult_type = rand(rng, [:Single, :Same, :Random])
    if mult_type == :Same && order > 1
        m_all = rand(rng, 2:order)
        m_cnt = div(multiples_cnt - 1, m_all) + 1
        multiples = ones(Int, m_cnt) * m_all
        multiples_cnt = length(multiples)
    elseif mult_type == :Random && order > 1
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
    else
        multiples = ones(Int, multiples_cnt)
    end
    uniques = sort(rand(rng, T, length(multiples)))
    axis = Glissa.axis_multiples_to_repeats(uniques, multiples)
    @assert length(axis) == sum(multiples)
    @assert length(uniques) == length(multiples)
    (axis, uniques, multiples, order)
end

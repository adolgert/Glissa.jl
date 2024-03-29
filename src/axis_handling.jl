using Random
using Distributions
using Logging


"""
    RepeatedIndex(multiple::AbstractVector)

Uses a multiplicity vector to present a vector as having repeated elements.
If a vector is length 3, and `multiple` is `[2, 1, 3]``, then this index
will return the values `[1, 1, 2, 3, 3, 3]`.
"""
struct RepeatedIndex{T<:Integer} <: AbstractUnitRange{T}
    # Store the cumulant instead of the multiplicity vector because it's
    # quicker to derive multiplicity from cumulant than to recompute the cumulant.
    cumulant::AbstractVector{Int}
    function RepeatedIndex{T}(multiple::AbstractVector) where {T<:Integer}
        return new(convert(Vector{Int}, cumsum(multiple)))
    end
end


function multiplicity(ri::RepeatedIndex, i)
    if i > 1
        ri.cumulant[i] - ri.cumulant[i-1]
    else
        ri.cumulant[1]
    end
end


function Base.iterate(ri::RepeatedIndex{T}) where {T}
    length(ri.cumulant) == 0 && return nothing
    iterate(ri, (1, 0))
end


function Base.iterate(ri::RepeatedIndex{T}, state) where T
    axis_idx, mult_idx = (state[1], state[2] + 1)
    while true
        mult_idx <= multiplicity(ri, axis_idx) && break
        axis_idx += 1
        axis_idx > length(ri.cumulant) && return nothing
        mult_idx = 1
    end
    axis_idx, (axis_idx, mult_idx)
end


# The a-th index of the multiples vector is from the i=\sum_{0}^{a-1} m_i entry,
# up to, but not including, the i=\sum_{0}^a m_i entry. This function inverts that,
# so, given a, find i.
function Base.getindex(ri::RepeatedIndex, i::Integer)
    searchsortedlast(ri.cumulant, i - 1) + 1
end


Base.first(::RepeatedIndex) = 1
Base.last(ri::RepeatedIndex) = ri.cumulant[end]
Base.IteratorSize(::Type{RepeatedIndex}) = HasLength()
Base.IteratorEltype(::Type{RepeatedIndex}) = HasEltype()
Base.eltype(::RepeatedIndex{T}) where {T} = T
Base.length(ri::RepeatedIndex) = length(ri.cumulant) > 0 ? ri.cumulant[end] : 0
# size(ri::RepeatedIndex) = (1,)


function axis_repeats_to_multiples(axis, ϵ = 0.0)
    uniques = copy(axis)
    multiplicity = zeros(Int, length(axis))
    ucnt = 0
    last = zero(eltype(axis))
    for uidx in eachindex(axis)
        if ucnt == 0 || abs(axis[uidx] - last) > ϵ
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
    ri = RepeatedIndex{Int}(multiplicity)
    axis = zeros(eltype(uniques), length(ri))
    for (aidx, midx) in enumerate(ri)
        axis[aidx] = uniques[midx]
    end
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
    sort(cover, dims=1, rev = false, lt=less_than_cover)
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

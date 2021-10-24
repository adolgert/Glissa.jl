function axis_repeats_to_multiples(axis)
    uniques = copy(axis)
    multiplicity = zeros(Int, length(axis))
    ucnt = 0
    last = axis[1] - one(eltype(axis))
    for uidx in 1:length(axis)
        if axis[uidx] != last
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
    @assert ax_idx = length(axis) + 1
    axis
end

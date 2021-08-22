"""
A recursive implementation of divided differences.
"""
mutable struct DividedDifference{T}
    memo::Dict{Tuple{Int64,Int64},T}
    DividedDifference{T}() where {T <: Real} = new(Dict{Tuple{Int64,Int64},T}())
end


"""
Calculate divided differences of a function f at points τ. This calculates the divided difference
from τ_i to τ_j. The `dd` is a memoization.
"""
function divided_difference(i::Integer, j::Integer, τ::Vector, f::Vector, dd::DividedDifference)
    if (i, j) ∈ keys(dd.memo)
        return dd.memo[(i, j)]
    end
    
    val = zero(Float64)
    if i == j
        val = f[i]
    elseif i < j
        val = (divided_difference(i + 1, j, τ, f, dd) - divided_difference(i, j - 1, τ, f, dd)) /
            (τ[j] - τ[i])
    else
        throw(ArgumentError("divided_difference needs i <= j, found i = {i}, j = {j}"))
    end
    dd.memo[(i, j)] = val
    val
end


function divided_difference_explicit(
    i::Integer, k::Integer, τ::Vector, f::Vector{T}) where {T <: Real}

    total = zero(T)
    for j = i:k
        denom = one(T)
        for l = i:k
            if l != j
                denom *= τ[j] -τ[l]
            end
        end
        total += f[j] / denom
    end
    total
end

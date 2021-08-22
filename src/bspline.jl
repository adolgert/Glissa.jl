
"""
An explicit representation of the normalized B-spline.
"""
function normalized_bspline(i::Integer, kp1::Integer, λ::Vector{T}, x::T) where {T <: Real}
    total = zero(T)
    for j = 0:kp1
        denom = one(T)
        for l = 0:(j - 1)
            denom *= λ[i + j] - λ[i + l]
        end
        for l = (j + 1):kp1
            denom *= λ[i + j] - λ[i + l]
        end
        total += truncated(λ[i + j], k, x) / denom
    end
    (λ[i + kp1] - λ[i]) * total
end

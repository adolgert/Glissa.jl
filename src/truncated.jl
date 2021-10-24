@doc raw"""
This is the truncated function that defines the Q B-spline.

$(x - c)_{+}^k$

It's only defined when $x > c$.
"""
function truncated(c, k, x)
    if x ≥ c
        (x - c)^k
    else
        zero(x)
    end
end


"""
A truncated function with a derivative with respect to `c`.
"""
function truncated(c, k, x, n)
    if x ≥ c
        if n == 0
            (x - c)^k
        elseif n ≤ k
            (-1)^n * prod((k-n+1):k) * (x - c)^(k - n)
        else
            zero(x)
        end
    else
        zero(x)
    end
end

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

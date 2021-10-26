# These algorithms are from Schumaker's book, Spline Functions: Basic Theory, 3rd ed.
# They generate the values, derivatives, and integrals of B-splines
# directly from an axis with its multiplicity. They don't create an intermediate
# representation of the polynomial for the B-spline.


@doc raw"""
    generate_normalized_bsplines!(N::AbstractVector{T}, y::AbstractVector, l, m, x::T) where {T}

Given x in ``[y_l, y_l+1)`` to generate ``N^m_{l+1-m}(x)`` to ``N_l^m(x).``
`y` is the axis. `m` is the order of the b-spline. The `l` is found by a search of the axis,
and the axis can have repeated elements. This algorithm applies to axis intervals where
$y_{l+1} > y_l$, not equal to it. That tells us how to number the expansion coefficients, `c`.
It's Algorithm 5.5 from Schumaker's book.
"""
function generate_normalized_bsplines!(
    N::AbstractVector{T}, y::AbstractVector, l, m, x::T) where {T}
    @assert length(N) == m + 1
    @assert length(y) > l

    Q = N  # Use the incoming storage as a buffer in which to calculate Q.
    Q[1:end] .= zero(T)
    if (y[l + 1] > y[l])
        Q[m] = one(T) / (y[l + 1] - y[l]) # else leave it zero.
    end
    for j in 2:(m - 1)
        # smallest i: m - (m-1) + 1, m -m, 1 + 1 = 2. largest = m
        for i in (m - j + 1):m
            denom = y[i + l - m + j] - y[i + l - m]
            a1 = (denom > 0) ? (x - y[i + l - m]) / denom : zero(T)
            a2 = one(T) - a1
            Q[i] = a1 * Q[i] + a2 * Q[i + 1]
        end
    end
    # This step converts b-splines Q into normalized b-splines N.
    for ii in 1:m
        N[ii] = (x - y[ii + l - m]) * Q[ii] + (y[ii + l] - x) * Q[ii + 1]
    end
    return nothing
end


# Algorithm 5.6: Evaluation of s(x) for Given a ≤ x < b.
# c are the expansion coefficients. m is the degrees of freedom, order + K.
# The expansion coefficients are numbered such that the first interval is 1-order.
# The next interval is 2-order, the next 3-order.
function evaluate_bspline56(c::AbstractArray{T}, y::AbstractArray, m, x::T) where {T}
    l = searchsortedlast(y, x)
    N = zeros(T, m + 1)
    generate_normalized_bsplines!(N, y, l, m, x)
    x = zero(T)
    for i = 1:m
        x += c[i + l - m] * N[i]
    end
    x
end


# Algorithm 5.8: Evaluate B-spline expansion at s Given x in [a,b]
function evaluate_bspline(c::AbstractArray{T}, y::AbstractArray, m, x::T) where {T}
    l = searchsortedlast(y, x)
    cx = similar(c)
    for init in 1:m
        cx[init] = c[init + l - m]
    end
    for j in 2:m
        for i in m:-1:j
            denom = y[i + l - j + 1] - y[i + l - m]
            a1 = (denom > 0) ? (x - y[i + l - m]) / denom : zero(T)
            a2 = one(T) - a1
            cx[i] = a1 * cx[i] + a2 * cx[i - 1]
        end
    end
    cx[m]
end


# Algorithm 5.10. Calculation of expansion coefficient matrix for derivatives.
# m is the order. c is the coefficient vector.
# Initialize `cd=zeros(T, m, length(c))`.
function derivative_expansion_coefficients!(cd, c::AbstractVector{T}, m) where T
    n = length(c)
    for init in 1:n
        cd[1, init] = c[init]
    end
    for j in 2:m
        mj = m - j + 1
        for i in n:-1:j
            denom = y[i + mj] - y[i]
            if denom ≠ 0
                cd[j, i] = mj * (cd[j - 1, i] - cd[j - 1, i - 1]) / denom
            else
                cd[j, i] = 0
            end
        end
    end
    nothing
end


# Algorithm 5.11 Compute D^{d-1}s(x) for given a ≤ x < b
# c is expansion coefficient. cd is derivatives matrix of expansion coefficient.
# m is order. d is derivative, x is value at which to evaluate.
function derivative_at(cd, m, d, x)
    l = searchsortedlast(y, x)
    return evaluate_bspline(cd[d, :], y, m - d + 1, x)
end


# Algorithm 5.12. Compute all derivatives
# Fill `s` with an array of all derivatives of the spline, starting from the 0th derivative.
# `cd` - the derivatives matrix of the expansion coefficient.
function all_derivatives!(s::AbstractVector, cd, m, x)
    l = searchsortedlast(y, x)
    s[m] = cd[m, l]
    for d = 1:(m-1)
        s[d] = evaluate_bspline(cd[d, :], y, m - d + 1, x)
    end
    nothing
end

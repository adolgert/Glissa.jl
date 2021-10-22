# These algorithms are from Schumaker's book, Spline Functions: Basic Theory, 3rd ed.
# They generate the values, derivatives, and integrals of B-splines
# directly from an axis with its multiplicity. They don't create an intermediate
# representation of the polynomial for the B-spline.

# Algorithm 5.3 5
# Given that x lies in [y_p, y_q) (That's closed on the left, open on the right)
# find l such that y_l <= x < y_{l+1}.
# This can return 0 if x < y[1]
function bisection(y::AbstractVector{T}, x::T) where {T}
    searchsortedlast(y, x)
end


# Algorithm 5.4 Given x in [y_p, y_q) and a guess for l,
# find l such that y_l <= x < y_l + 1.
function find_left_index(y::AbstractVector{T}, x::T, l::Int) where {T}
    n = length(y)
    if l + 1 ≤ n && x ≥ y[l + 1]
        if l + 1 == n
            return l + 1
        elseif l + 2 ≤ n && x < y[l + 2]
            return l + 1
        elseif l + 3 ≤ n && x < y[l + 3]
            return l + 2
        else
            return l + 2 + searchsortedlast(view(y, (l + 3):n), x)
        end
    elseif x < y[l]
        if l - 1 > 0
            if y[l - 1] ≤ x
                return l - 1
            else
                return searchsortedlast(view(y, 1:(l - 1)), x)
            end
        else
            return 0
        end
    else
        return l
    end
end


# Algorithm 5.5: Given x in [y_l, y_l+1) to generate N^m_{l+1-m}(x) to N_l^m(x).
# y is the axis. m is the order of the b-spline. The $l$ is found by a search of the axis,
# and the axis can have repeated elements. This algorithm applies to axis intervals where
# $y_{l+1} > y_l$, not equal to it. That tells us how to number the expansion coefficients, $c$.
function generate_normalized_bsplines!(
    N::AbstractVector{T}, y::AbstractVector, l, m, x::T) where {T}
    @assert length(N) == m + 1
    @assert length(y) > l

    Q = N  # Use the incoming storage as a buffer in which to calculate Q.
    Q[1:(m-1)] .= zero(T)
    Q[m] = (y[l + 1] > y[l]) ? one(T) / (y[l + 1] - y[l]) : zero(T)
    Q[m + 1] = zero(T)
    for j in 2:(m - 1)
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


@doc raw"""
Algorithm 5.15: Construct a piecewise polynomial representation of a polynomial spline
that is represented by spline coefficients. This converts out of B-splines to the polynomial
representation. If there are more than about 2 evaluations of each spline interval, this
can be more efficient.
Each spline ``i`` is represented by a polynomial of the form

``s_i(x)= \sum_{j=1}^m w_{ji} (x - \tau_i)^{j-1}``

where ``w_{ij}`` are polynomial constants for the `i`th B-spline.

Initialize `w = zeros(T, m, n)` where `n` is the number of B-spline coefficients.
"""
function piecewise_representation!(
      w::AbstractArray{T}, cd::AbstractArray{T}, τ::AbstractVector) where {T}

    m, n = size(cd)  # m is the order. n is the number of B-splines.
    s = zeros(T, m)
    factorials = cumprod(Iterators.flatten((1:1, 1:(m-1)))) # factorial(j-1)
    for i = 2:n
        all_derivatives!(s, cd, m, τ[i])
        for j = 1:m
            w[j, i] = s[j] / factorials[j]
        end
    end
    # this is for x_1 from below.
    all_derivatives!(s, cd, m, τ[1])
    for jj = 1:m
        w[jj, 1] = s[jj] / factorial(j - 1)
    end
    nothing
end

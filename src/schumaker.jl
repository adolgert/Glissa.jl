# These algorithms are from Schumaker's book, Spline Functions: Basic Theory, 3rd ed.

# Algorithm 5.3 Bisection
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
# y is the axis. m is the order of the b-spline.
function generate_normalized_bsplines!(
    N::AbstractVector{T}, y::AbstractVector{T}, l, m, x) where {T}
    @assert length(N) == m + 1
    @assert length(y) > l

    Q = N
    Q[m] = (y[l + 1] > y[l]) ? 1 / (y[l + 1] - y[l]) : 0
    Q[m + 1] = 0
    for j in 2:(m - 1)
        for i in (m - j + 1):m
            denom = y[i + l - m + j] - y[i + l - m]
            a1 = (denom > 0) ? (x - y[i + l - m]) / denom : 0
            a2 = 1 - a1
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
function evaluate_bspline56(c, y, m, x)
    l = searchsortedlast(y, x)
    N = zeros(T, m + 1)
    generate_normalized_bsplines!(N, y, l, m, x)
    x::T = 0
    for i = 1:m
        x += c[i] * N[i]
    end
    x
end


# Algorithm 5.8: Evaluate B-spline expansion at s Given x in [a,b]
function evaluate_bspline(c, y, m, x)
    l = searchsortedlast(y, x)
    cx = similar(c)
    for init in 1:m
        cx[init] = c[init + l - m]
    end
    for j in 2:m
        for i in m:-1:j
            denom = y[i + l - j + 1] - y[i + l - m]
            a1 = (x - y[i + l - m]) / denom
            a2 = 1 - a1
            cx[i] = a1 * cx[i] + a2 * cx[i - 1]
        end
    end
    cx[m]
end


# Algorithm 5.10. Calculation of expansion coefficient matrix for derivatives.
# m is the order. c is the coefficient vector.
function derivative_expansion_coefficients(c::AbstractVector{T}, m) where T
    n = length(c)
    cd = zeros(T, m, n)
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
    cd
end


# Algorithm 5.11 Compute D^{d-1}s(x) for given a ≤ x < b
# c is expansion coefficient. cd is derivatives matrix of expansion coefficient.
# m is order. d is derivative, x is value at which to evaluate.
function derivative_at(cd, m, d, x)
    l = searchsortedlast(y, x)
    return evaluate_bspline(cd[d, :], y, m - d + 1, x)
end


# Algorithm 5.12. Skipping for now. Computes all derivatives.

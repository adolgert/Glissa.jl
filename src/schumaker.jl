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

@doc raw"""
This is a piecewise polynomial. It has N+1 points on the abcissa, ``\tau``.
It has N cubic interpolants in `c`, which is size 4 x N for a cubic spline.
At each interval, the polynomial is in Horner Form,

``P_i(x) = c_{1,i} + c_{2,i}(x-x_i) + c_{3,i}(x-x_i)^2 + c_{4,i}(x-x_i)^3``

This spline has one type for the abcissa and one for the ordinate coefficients.
It should be the case that one(T) * one(X) is of type T.

This can be a spline of another order. It's just a matter of the dimensions
of the coefficient array, `c`. You would evaluate this with
`evalpoly(x-τ[i], c[1:end, i])`, where `i` is the interval.

Every piecewise polynomial is a functor, so you can evaluate it
at a point with `ps(x)`.
"""
abstract type PiecewisePolynomial{X,T} end


struct PolySpline{X,T} <: PiecewisePolynomial{X,T}
    τ::AbstractVector{X}  # The abcissa
    c::AbstractMatrix{T}  # Polynomial coefficients (power of x, interval index).
end

"""
    order(ps::PiecewisePolynomial)

The order is one plus the degree of the polynomial.
"""
order(ps::PiecewisePolynomial) = size(ps.c, 1)

@doc raw"""
    degree(ps::PiecewisePolynomial)

Polynomial degree, so 2 for ``x^2``
"""
degree(ps::PiecewisePolynomial) = order(ps) - 1

"""
    length(ps::PiecewisePolynomial)

Number of vertices in axis, which is 1 + number of intervals
"""
Base.length(ps::PiecewisePolynomial) = length(ps.τ)

"""
    eltype(::PiecewisePolynomial{X,T})

Type of the polynomial constants.
"""
Base.eltype(::PiecewisePolynomial{X,T}) where {X,T} = T


# Difficult to add this document string with Documenter.jl.
# Evaluate a piecewise polynomial at `x`.
# """
function (cs::PiecewisePolynomial)(x::A) where {A <: Real}
    # Index of first greater-than-or-equal-to x.
    i = Sort.searchsortedlast(cs.τ, x)
    # For points off the sides, use the closest polynomial.
    i = max(min(length(cs.τ) - 1, i), 1)
    evalpoly(x - cs.τ[i], view(cs.c, :, i))
end

"""
    evaluate!(cs::PiecewisePolynomial, x::AbstractVector, y::AbstractVector)

Evaluate a piecewise polynomial at a series of `x` values and write into `y`.
"""
function evaluate!(cs::PiecewisePolynomial, x::AbstractVector, y::AbstractVector)
    # Index of first greater-than-or-equal-to x.
    k = 1
    i = Sort.searchsortedfirst(cs.τ, x[k]) - 1
    i = max(min(length(cs.τ) - 1, i), 1)
    while k <= length(x)
        y[k] = evalpoly(x[k] - cs.τ[i], view(cs.c, :, i))
        k += 1
        while i + 1 < length(cs.τ) && cs.τ[i + 1] < x[k]
            i += 1
        end
    end
end

"""
    integral_in_interval(cs::PiecewisePolynomial{X,T}, i, x::A) where {A <: Real, X, T}

Integrate a spline from a knot at `i` to the value `x`.
# (((A(1/4)*cs.c[4, i] * Δ + A(1/3)*cs.c[3, i]) * Δ + A(1/2)*cs.c[2, i]) * Δ + cs.c[1, i])*Δ
"""
function integral_in_interval(cs::PiecewisePolynomial{X,T}, i, x::A) where {A <: Real, X, T}
    Δ = x - cs.τ[i]
    m = size(cs.c, 1)
    if m > 0
        total = zero(T)
        for d = m:-1:1
            total += cs.c[d, i] / T(d)
            total *= Δ
        end
        total
    else
        zero(T)
    end
end


"""
    integrate(cs::PiecewisePolynomial, x1::A, x2::A) where {A <: Real}

Integrate a spline from `x1` to `x2`.
"""
function integrate(cs::PiecewisePolynomial, x1::A, x2::A) where {A <: Real}
    # Index of first greater-than-or-equal-to x.
    i = Sort.searchsortedfirst(cs.τ, x1) - 1
    # For points off the sides, use the closest polynomial.
    i = max(min(length(cs.τ) - 1, i), 1)
    j = Sort.searchsortedfirst(cs.τ, x2) - 1
    # For points off the sides, use the closest polynomial.
    j = max(min(length(cs.τ) - 1, j), 1)
    intervening_intervals = zero(A)
    for k in i:(j-1)
        intervening_intervals += integral_in_interval(cs, k, cs.τ[k + 1])
    end
    intervening_intervals + integral_in_interval(cs, j, x2) - integral_in_interval(cs, i, x1)
end


"""
    derivative(cs::PiecewisePolynomial{X,T}) where {X,T}

Given a spline, create a new spline that is its derivative.
"""
function derivative(cs::PiecewisePolynomial{X,T}) where {X,T}
    m = size(cs.c, 1)
    c2 = similar(cs.c[2:m, :])
    for i = 1:size(cs.c, 2)
        for d = 1:(m - 1)
            c2[d, i] = T(d) * cs.c[d + 1, i]
        end
    end
    PolySpline{X,T}(cs.τ, c2)
end


"""
    Base.:*(cs1::PiecewisePolynomial{X,T}, cs2::PiecewisePolynomial{X,T}) where {X,T}

Multiply two splines of order `m1` and `m2` to get a spline of order `m1+m2`.
This only works if the two splines have the same abcissa.
"""
function Base.:*(cs1::PiecewisePolynomial{X,T}, cs2::PiecewisePolynomial{X,T}) where {X,T}
    m1 = size(cs1.c, 1)
    m2 = size(cs2.c, 1)
    c3 = zeros(T, m1 + m2 - 1, size(cs1.c, 2))
    for i in 1:size(cs1.c, 2)
        for d1 = 1:m1
            for d2 = 1:m2
                c3[d1 + d2 - 1, i] += cs1.c[d1, i] * cs2.c[d2, i]
            end
        end
    end
    PolySpline{X, T}(cs1.τ, c3)
end

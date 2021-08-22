using Base.Sort
using LinearAlgebra

@doc raw"""
This is a traditional cubic spline. It has N+1 points on the abcissa, ``\tau``.
It has N cubic interpolants in `c`, which is size 4 x N.
At each interval, the polynomial is in Horner Form,

``P_i(x) = c_{1,i} + c_{2,i}(x-x_i) + c_{3,i}(x-x_i)^2 + c_{4,i}(x-x_i)^3``

This spline has one type for the abcissa and one for the ordinate coefficients.
It should be the case that one(T) * one(X) is of type T.
"""
struct CubicSpline{X,T}
    τ::Vector{X}
    c::Matrix{T}
end


function horner_in_interval(cs::CubicSpline, i, x)
    Δ = x - cs.τ[i]
    ((cs.c[4, i] * Δ + cs.c[3, i]) * Δ + cs.c[2, i]) * Δ + cs.c[1, i]
end


function (cs::CubicSpline)(x::A) where {A <: Real}
    # Index of first greater-than-or-equal-to x.
    i = Sort.searchsortedfirst(cs.τ, x) - 1
    # For points off the sides, use the closest polynomial.
    i = max(min(length(cs.τ) - 1, i), 1)
    horner_in_interval(cs, i, x)
end


function evaluate!(cs::CubicSpline, x::Vector, y::Vector)
    # Index of first greater-than-or-equal-to x.
    k = 1
    i = Sort.searchsortedfirst(cs.τ, x[k]) - 1
    i = max(min(length(cs.τ) - 1, i), 1)
    while k <= length(x)
        y[k] = horner_in_interval(cs, i, x[k])
        k += 1
        while i + 1 < length(cs.τ) && cs.τ[i + 1] < x[k]
            i += 1
        end
    end
end


@doc raw"""
τ are the abcissa from 1 to N+1. f are the known values of length N+1.
fp_bounds are derivatives at the first and last endpoint.
fp is the output, of length N + 1.
This is from Conte and deBoor Sec. 4.8, Eq. 4.60. The equation is this, with ``i``
running from 2 to N.

``(\Delta x_i) s_{i-1}+2(\Delta x_{i-1}+\Delta x_i)s_i +(\Delta x_{i-1}s_{i+1}) =
3(f[x_{i-1},x_i]\Delta x_i+f[x_i,x_{i+1}]\Delta x_{i-1})``

Here, ``\Delta x_i = x_{i+1} -x{i}`` and

``f[x_i,x_{i+1}] = \frac{f[x_{i+1}] - f[x_i]}{x_{i+1}-x_i}``

"""
function global_derivatives!(τ, f::Vector{T}, fp) where {T <: Real}
    # This implementation builds the tridiagonal matrix and then solves it.
    # You could write a Gaussian elimination to generate terms on the fly, saving memory.
    # There are N-1 equations.
    N = length(τ) - 1
    # lower, mid, and upper of tridiagonal matrix.
    dl = zeros(T, N - 2)  # indices from 3 to N
    d = zeros(T, N - 1) # indices from 2 to N
    du = zeros(T, N - 2) # indices from 2 to N-1
    rhs = zeros(T, N - 1) # indices from 2 to N

    for il = 3:N
        dl[il - 2] = τ[il + 1] - τ[il]  # \Delta x_i
    end
    for im = 2:N
        d[im - 1] = T(2) * (τ[im + 1] - τ[im - 1])  # 2(\Delta x_{i-1} + \Delta x_i)
    end
    for il = 2:(N - 1)
        du[il - 1] = τ[il] - τ[il - 1]  # \Delta x_{i-1}
    end
    for il = 2:N
        fleft = (f[il] - f[il - 1]) / (τ[il] - τ[il - 1])  # f[x_{i-1},x_i]
        fright = (f[il + 1] - f[il]) / (τ[il + 1] - τ[il])  # f[x_i,x_{i+1}]
        # 3(f[x_{i-1},x_i]\Delta x_i+f[x_i,x_{i+1}]\Delta x_{i-1})
        rhs[il - 1] = T(3) * (fleft * (τ[il + 1] - τ[il]) + fright * (τ[il] - τ[il - 1]))
    end
    # The known derivatives at the boundaries get moved from the left-hand side to the RHS.
    rhs[1] -= fp[1] * (τ[3] - τ[2])  # i = 2, (\Delta x_i) s_{i-1}
    rhs[N - 1] -= fp[end] * (τ[N] - τ[N - 1])  # i = N, \Delta x_{i-1} s_{i+1}
    tmat = Tridiagonal(dl, d, du)
    fp[2:(end - 1)] .= tmat \ rhs
    nothing
end


"""
Given an abcissa, τ, of length N+1 values f of length N+1, and derivatives s
of length N+1, this computes the 4xN matrix c of cubic coefficients.
Conte and deBoor Eq. 4.55
"""
function cubic_spline_coefficients!(τ, f::Vector{T}, s, c) where {T <: Real}
    for i = 1:(length(τ) - 1)
        c[1, i] = f[i]
        c[2, i] = s[i]
        Δ = τ[i + 1] - τ[i]
        dd = (f[i + 1] - f[i]) / Δ
        c4p = (s[i + 1] + s[i] - T(2) * dd) / Δ
        c[3, i] = (dd - s[i]) / Δ - c4p
        c[4, i] = c4p / Δ
    end
    nothing
end


function cubic_spline_flat_endpoints(τ::Vector{X}, f::Vector{T}) where {X <: Real, T <: Real}
    N = length(τ) - 1
    s = zeros(T, N + 1)
    s[1] = zero(T)
    s[end] = zero(T)
    global_derivatives!(τ, f, s)
    c = zeros(T, 4, N)
    cubic_spline_coefficients!(τ, f, s, c)
    CubicSpline{X,T}(τ, c)
end


"""
"Accurate Monotonicity Preserving Cubic Interpolation" by James M Hyman. 1983.
"""
function deboor_swartz_criterion(s::T, sm1, sp1) where {T <: Real}
    smin = min(sm1, sp1)
    smax = max(sm1, sp1)
    # Equation 2.3
    if zero(T) < smin
        s = min(max(zero(T), s), T(3) * smin)
    elseif smax < zero(T)
        s = max(min(zero(T), s), T(3) * smax)
    elseif sm1 * sp1 <= zero(T)
        s = zero(T)
    else
        s  # Data are locally monotone and within bounds.
    end
end

function hyman_criterion(s::T, sm1, sp1) where {T <: Real}
    # Equation 2.6 "extends" equation 2.3.
    smin = min(sm1, sp1)
    smax = max(sm1, sp1)
    if sm1 * sp1 > 0
        s = sp1
    end
    if s >= 0
        s = min(max(0, s), T(3) * min(abs(smin), abs(smax)))
    else
        s = max(min(0, s), -T(3) * min(abs(smin), abs(smax)))
    end
    s
end


"""
Uses Hyman piecewise monotonicity.
"""
function project_to_monotonicity!(x, f, s)
    s[1] = hyman_criterion(s[1], (f[2] - f[1])/(x[2] - x[1]), (f[2] - f[1])/(x[2] - x[1]))
    for i in 2:(length(s)-1)
        s[i] = hyman_criterion(
            s[i],
            (f[i+1] - f[i])/(x[i+1] - x[i]),
            (f[i] - f[i-1])/(x[i] - x[i-1])
            )
    end
    s[end] = hyman_criterion(
        s[end],
        (f[end] - f[end-1])/(x[end]-x[end-1]),
        (f[end] - f[end-1])/(x[end]-x[end-1])
        )
    nothing
end

# Options for endpoints
struct ZeroDerivativeEndpoints
end

struct FlatEndpoints
end

# Options for monotonicity.
struct FreeFit
end

struct MonotonicFit
end

# Options for determination of derivatives
struct GlobalDerivatives
end

function cubic_spline(τ::Vector{X}, f::Vector{T}, fp::Vector{T}) where {X <: Real, T <: Real}
    N = length(τ) - 1
    s = zeros(T, N + 1)
    s[1] = fp[1]
    s[end] = fp[2]
    global_derivatives!(τ, f, s)
    project_to_monotonicity!(τ, f, s)
    c = zeros(T, 4, N)
    cubic_spline_coefficients!(τ, f, s, c)
    CubicSpline{X,T}(τ, c)
end

using Base.Sort
using LinearAlgebra

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
function global_derivatives!(τ, f::AbstractVector{T}, fp) where {T <: Real}
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
    fpans = tmat \ rhs
    fp[2:(end - 1)] .= fpans
    nothing
end


"""
Given an abcissa, τ, of length N+1 values f of length N+1, and derivatives s
of length N+1, this computes the 4xN matrix c of cubic coefficients.
Conte and deBoor Eq. 4.55
"""
function cubic_spline_coefficients!(τ, f::AbstractVector{T}, s, c) where {T <: Real}
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


"""
This comes from Hyman's paper but is modified by the R implementation of splinefun,
which Simon Wood wrote. I didn't understand from the paper the exact structure of
decisions to take in the function.
"""
function hyman_criterion(s::T, sm1, sp1) where {T <: Real}
    # Equation 2.6 "extends" equation 2.3.
    # The paper doesn't specify setting s=sp1 _before_ the if-then about the sign of s.
    σ = s
    if sm1 * sp1 > 0
        σ = sp1  # The paper seems to say this should be set to zero.
    end
    if σ >= 0
        s = min(max(0, s), T(3) * min(abs(sm1), abs(sp1)))
    else
        s = max(min(0, s), -T(3) * min(abs(sm1), abs(sp1)))
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
            (f[i] - f[i-1])/(x[i] - x[i-1]),
            (f[i+1] - f[i])/(x[i+1] - x[i])
            )
    end
    s[end] = hyman_criterion(
        s[end],
        (f[end] - f[end-1])/(x[end]-x[end-1]),
        (f[end] - f[end-1])/(x[end]-x[end-1])
        )
    nothing
end


"""
This is how splinefun sets the coefficients. Seems to have the same
result as the code above.
"""
function hyman_coefficients!(τ, f::AbstractVector{T}, s, c) where {T <: Real}
    for i = 1:(length(τ) - 1)
        c[1, i] = f[i]
        c[2, i] = s[i]
        if i < length(τ)
            Δ = τ[i + 1] - τ[i]
            y = -(f[i + 1] - f[i])
            c[3, i] = -(3y + (2s[i] + s[i+1])Δ) / Δ^2
            c[4, i] = (2y/Δ + s[i] + s[i+1])/Δ^2
        else
            # That's correct. This doesn't get called, but R's splinefun sets these.
            # Leaving here until I figure out why we need a polynomial defined for the
            # region to the right of the last input point.
            Δ = τ[i] - τ[i-1]
            y = -(f[i] - f[i-1])
            c[3, i] = (3y + (s[i-1] + 2s[i])Δ) / Δ^2
            c[4, i] = c[4, i-1]
        end
    end
    nothing
end


# Options for endpoints
struct ZeroDerivativeEndpoints
end

struct FlatEndpoints
end

# Options for monotonicity.
struct FreeSlope
end

function project_slope!(::FreeSlope, τ, f, s)
end


struct Monotonic
end

function project_slope!(::Monotonic, τ, f, s)
    project_to_monotonicity!(τ, f, s)
end


# Options for determination of derivatives
struct GlobalDerivatives
end


function cubic_spline(
    τ::AbstractVector{X}, f::AbstractVector{T}, fp; slope=FreeSlope()
    ) where {X <: Real, T <: Real}

    cubic_spline(τ, f, convert(typeof(f), fp); slope = slope)
end


function cubic_spline(
    τ::AbstractVector{X}, f::AbstractVector{T}, fp::AbstractVector{T};
    slope=FreeSlope()
    ) where {X <: Real, T <: Real}

    N = length(τ) - 1
    s = zeros(T, N + 1)
    s[1] = fp[1]
    s[end] = fp[2]
    global_derivatives!(τ, f, s)
    project_slope!(slope, τ, f, s)
    c = zeros(T, 4, N)
    cubic_spline_coefficients!(τ, f, s, c)
    PolynomialSpline{X,T}(τ, c)
end

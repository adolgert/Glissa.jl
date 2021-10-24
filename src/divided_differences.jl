using LinearAlgebra

# This file is about calculating divided differences.
# Divided differences are the most common first presentation of B-splines.
# When the axis has no repeated knots, these are simple to calculate, but I'm
# having trouble finding clear descriptions of divided differences when the
# knots repeat.
#
# You will see below several versions of definitions, so that I can check them
# against each other.
#
# 1. Ratio of determinants for unique axis knots.
# 2. Ratio of determinants for repeated axis knots.
# 3. Recursive definition, applies to repeated knots.
# 4. Explicit formula, only for unique knots.
# 5. The B-splines themselves, defined using these divided differences.


# Build divided-differences from this matrix representation.
# This is a very theoretical way to build a divided-difference, which makes it less
# susceptible to accidental typos.
# u1(t1) u2(t1) u3(t1)  t are rows. u are columns.
# u1(t2) u2(t2) u3(t2)
# u1(t3) u2(t3) u3(t3)
# The cross-matrix represents a matrix of functions over the axis.
# This makes an AbstractMatrix so that we can call the usual LinearAlgebra.det to
# get the determinant. Just makes it easier.
struct CrossMatrix{T} <: AbstractMatrix{T}
    t::AbstractVector  # An array of axis coordinates
    u::AbstractVector  # An array of functions
end

function cross_matrix(axis, ufuncs)
    T = typeof(ufuncs[1](axis[1]))
    CrossMatrix{T}(axis, ufuncs)
end

Base.getindex(cm::CrossMatrix, i, j) = cm.u[j](cm.t[i])
Base.size(cm::CrossMatrix) = (length(cm.t), length(cm.u))


function vandermonde_determinant(axis)
    total = zero(eltype(axis))
    for i in 1:(length(axis) - 1)
        for j in (i+1):length(axis)
            total *= axis[j] - axis[i]
        end
    end
    total
end


@doc raw"""The equation 2.67 in Schumaker is hard to understand. It reads, exactly,

$V(t_1,\ldots,t_m)=\prod_{1\le i<j\le d}(\tau_j-\tau_i)^{l_jl_i} \prod_{i=1}^d\prod_{\nu=1}^{l_i-1}\nu!$

The LHS goes to $m$ to count the total length of the axis. The RHS uses $d$ which is the number
of unique values. I have a hard time believing it's a product of factorials at the end.
Maybe he wanted just a single product over the factorial?

$\prod_{i=1}^d (\l_i-1)!$

Let's check it.
"""
function vandermonde_determinant_repeats(uniques, multiplicity)
    total = zero(eltype(uniques))
    ucnt = length(uniques)

    for i in 1:(ucnt - 1)
        for j in (i+1):ucnt
            total *= (uniques[j] - uniques[i])^(multiplicity[i] * multiplicity[j])
        end
    end
    for fidx in 1:ucnt
        for midx in 1:(multiplicity[fidx] - 1)
            total *= factorial(midx)
        end
    end
    total
end


@doc raw"""
Divided differences that allow for derivatives.
Schumaker page 46-7, eqn 2.89.
"""
function divdiff_repeats(axis, f)
    uniques, multiples = axis_repeats_to_multiples(axis)
    total = zero(eltype(axis))
    denominator = vandermonde_determinant_repeats(uniques, multiplicity)
    for ax_idx in 1:length(uniques)
        msave = multiplicity[ax_idx]
        for der_idx in 1:msave
            multiplicity[ax_idx] = der_idx - 1  # This can be zero.
            numerator = vandermonde_determinant_repeats(uniques, multiplicity)
            if der_idx == 1
                total += numerator * f(uniques[ax_idx])
            else
                # The extra argument gets you a derivative.
                total += numerator * f(uniques[ax_idx], der_idx - 1)
            end
        end
        multiplicity[ax_idx] = msave
    end
    total / denominator
end


@doc raw"""
Definition 2.49 of divided differences from Schumaker. Eqn 2.86.

> axis = 0:5
> f(x) = 2.3 + 0.5*x^2 - x^3
> divdiff(axis, f)
"""
function divdiff(axis, f)
    u = vcat([x->x^(j-1) for j in 1:(length(axis) - 1)], [f])
    numerator = cross_matrix(axis, u)
    denominator = cross_matrix(axis, [x->x^(j-1) for j in 1:length(axis)])
    LinearAlgebra.det(numerator) / LinearAlgebra.det(denominator)
end


"""
A recursive implementation of divided differences.
"""
mutable struct DividedDifference{T}
    memo::Dict{Tuple{Int64,Int64},T}
    DividedDifference{T}() where {T <: Real} = new(Dict{Tuple{Int64,Int64},T}())
end


# Compute the derivative of a function with respect to the axis.
function divdiff_der(::Function, n)
    if n > 0
        error("no derivative defined")
    else
        f
    end
end


@doc raw"""
Calculate divided differences of a function f at points τ, using recursion.

$[\tau_i, \tau_{i+1},\ldots,\tau_j]f$

This calculates the divided difference from τ_i to τ_j. The `dd` is a memoization.
If any knots are equal, define `f(x,n)` where `n` is how many times to take the derivative.
"""
function divided_difference(i::Integer, j::Integer, τ::AbstractVector, f::Function, dd::DividedDifference)
    if (i, j) ∈ keys(dd.memo)
        return dd.memo[(i, j)]
    end
    
    val = zero(eltype(τ))
    if i == j
        val = f(τ[i])
    elseif i < j
        if all(τ[i:j] .== τ[i])
            val = f(τ[i], nth_derivative) / factorial(nth_derivative)
        else
            val = (divided_difference(i + 1, j, τ, f, dd) - divided_difference(i, j - 1, τ, f, dd)) /
                (τ[j] - τ[i])
        end
    else
        throw(ArgumentError("divided_difference needs i <= j, found i = {i}, j = {j}"))
    end
    dd.memo[(i, j)] = val
    val
end


"""
This explicit divided difference applies only when the axis is distinct.
"""
function divided_difference_explicit(
    i::Integer, k::Integer, τ::AbstractVector, f::Function)

    for check_idx in i:(k-1)
        @assert τ[check_idx] < τ[check_idx + 1]
    end

    T = typeof(f(τ[i]))
    total = zero(T)
    for j = i:k
        denom = one(T)
        for l = i:k
            if l != j
                denom *= τ[j] -τ[l]
            end
        end
        total += f(τ[j]) / denom
    end
    total
end


@doc raw"""
A Q version of the B-spline on an axis y, order m, index i, at x.
This version is defined using divided differences.
"""
function splineq_divided(y, m, i, x)
    if y[i] ≤ x < y[i + m]
        f(y, n = 0) = truncated(y, m - 1, x, n)
        memo = DividedDifference{eltype(y)}()
        (-1)^m * divided_difference(i, i + m, y, f, memo)
    else
        zero(eltype(y))
    end
end


@doc raw"""
An N version of the B-spline on an axis y, order m, index i, at x.
This version is defined using divided differences.
"""
function splinen_divided(y, m, i, x)
    (y[i + m] - y[i]) * splineq_divided(y, m, i, x)
end

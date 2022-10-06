# Reference

```@meta
CurrentModule = Glissa
```

## Axis Handling

These are helper functions to convert between different representations of
the ordinal axis because different authors refer to knots in different ways.

```@docs
RepeatedIndex
iterate(ri::RepeatedIndex{T}) where {T}
axis_repeats_to_multiples(axis, ϵ = 0.0)
axis_multiples_to_repeats(uniques, multiplicity)
bspline_count(axis::AbstractVector{T}, order) where {T <: Real}
MultIndex
bspline_indices_in_interval(multiples, order)
generate_random_polyspline_specification(rng, T)
```

## Piecewise Polynomial Representation

This is the most basic representation of a piecewise polynomial. This
module will use it to test spline solutions.

```@docs
PiecewisePolynomial
order(ps::PiecewisePolynomial)
degree(ps::PiecewisePolynomial)
Base.length(ps::PiecewisePolynomial)
Base.eltype(::PiecewisePolynomial)
evaluate!(cs::PiecewisePolynomial, x::AbstractVector, y::AbstractVector)
integral_in_interval(cs::PiecewisePolynomial{X,T}, i, x::A) where {A <: Real, X, T}
integrate(cs::PiecewisePolynomial, x1::A, x2::A) where {A <: Real}
derivative(cs::PiecewisePolynomial{X,T}) where {X,T}
Base.:*(cs1::PiecewisePolynomial{X,T}, cs2::PiecewisePolynomial{X,T}) where {X,T}
```

## Newton Form

The [Newton Polynomial](https://en.wikipedia.org/wiki/Newton_polynomial) represents
a polynomial on an interval by the values it takes on that interval. It's a particularly
simple expression of polynomials.

```@docs
NewtonForm
NewtonForm(τ::Vector, f::Vector{T}) where {T <: Real}
single_polynomial(τ::Vector{X}, f::Vector{T}) where {T <: Real, X<: Real}
```

## Truncated Linear Functions

B-splines are most commonly defined in terms of
[divided differences](https://en.wikipedia.org/wiki/Divided_differences). Those divided
differences define the coefficients with which to multiply the truncated functions
defined in this section of the library.

```@docs
truncated(c, k, x)
truncated(c, k, x, n)
normalized_bspline(i::Integer, kp1::Integer, λ::Vector{T}, x::T) where {T <: Real}
```

## Divided Differences

This section is about calculating divided differences.
Divided differences are the most common first presentation of B-splines.
When the axis has no repeated knots, these are simple to calculate, but I'm
having trouble finding clear descriptions of divided differences when the
knots repeat.

You will see below several versions of definitions, so that I can check them
against each other.

1. Ratio of determinants for unique axis knots.
2. Ratio of determinants for repeated axis knots.
3. Recursive definition, applies to repeated knots.
4. Explicit formula, only for unique knots.
5. The B-splines themselves, defined using these divided differences.

```@docs
CrossMatrix
vandermonde_determinant(axis)
vandermonde_determinant_repeats(uniques, multiplicity)
divdiff_repeats(axis, f)
divdiff(axis, f)
DividedDifference{T}
divided_difference(i::Integer, j::Integer, τ::AbstractVector, f::Function, dd::DividedDifference)
divided_difference_explicit(i::Integer, k::Integer, τ::AbstractVector, f::Function)
splineq_divided(y, m, i, x)
splinen_divided(y, m, i, x)
```

## Recursive B-splines

B-splines were defined recursively first.

```@docs
bsplineq_recursive(y, m, i, x)
bsplineq_recursive(y, m, i, x, n)
```

## Solving for B-Splines by Constraining Piecewise Splines


These functions solve for B-splines of order 2 by constructing a matrix
of constraints and solving that matrix for coefficients of each part of the
polynomial. This will give you coefficients for any B-spline. If you want
to know coefficients for B-splines on an integer axis of a given order, then
run this with the Rational type, and it should agree with published examples.

This isn't used, in practice, when working with B-splines, because you don't need
to represent them as polynomials. Instead, there are other functions that evaluate
the B-spline directly from the axis and smoothness conditions that define it.

```@docs
bspline_by_matrix!
reduce_axis(multiplicity::AbstractVector, order, i)
polyspline_constraints!(axis::AbstractArray{T}, multiplicity::AbstractVector{Int}, order) where {T}
```

## B-Splines (Official version)

These algorithms are from Schumaker's book, Spline Functions: Basic Theory, 3rd ed.
They generate the values, derivatives, and integrals of B-splines
directly from an axis with its multiplicity. They don't create an intermediate
representation of the polynomial for the B-spline.

```@docs
generate_normalized_bsplines!
evaluate_bspline56(c::AbstractArray{T}, y::AbstractArray, m, x::T) where {T}
evaluate_bspline(c::AbstractArray{T}, y::AbstractArray, m, x::T) where {T}
derivative_expansion_coefficients!(cd, c::AbstractVector{T}, m) where T
derivative_at(cd, m, d, x)
all_derivatives!(s::AbstractVector, cd, m, x)
```

## Cubic Splines

This file is a traditional cubic splines interpolation, including
monotonic splines.

```@docs
global_derivatives!(τ, f::AbstractVector{T}, fp) where {T <: Real}
cubic_spline_coefficients!(τ, f::AbstractVector{T}, s, c) where {T <: Real}
deboor_swartz_criterion(s::T, sm1, sp1) where {T <: Real}
hyman_criterion(s::T, sm1, sp1) where {T <: Real}
project_to_monotonicity!(x, f, s)
hyman_coefficients!(τ, f::AbstractVector{T}, s, c) where {T <: Real}
ZeroDerivativeEndpoints
FlatEndpoints
FreeSlope
Monotonic
cubic_spline
```

## Schumaker

These are other algorithms from Schumaker's book.

```@docs
piecewise_representation!
```

## Index

```@index
```

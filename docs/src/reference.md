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

## Index

```@index
```

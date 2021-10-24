# Test Plan

This package provides non-uniform and uniform B-spline implementations. The standard methods, from Schumaker, Dierkx, and back to Conte and de Boor, are complicated enough, but we might implement signal processing techniques, as well. All need to be checked.

## Risks for Testing

This isn't a list of project risks but of where I see complexity in the code.

1. Differences in notation among source material mean I hook together functions incorrectly.
1. Off-by-one errors coming from translation among C, Fortran, and Julia.
1. There could be domain problems when knots are very close together. I don't know how well-tested this published code is. Often, professional solutions don't look like the published ones, and the published ones lead to floating-point exceptions.
1. There could be interactions among types in Julia, that aren't anticipated. For instance, using Rationals with Floats, or Floats of differing precision.
1. Current usage doesn't exercise the code enough to see errors, for instance the use of non-uniform axes that also have inhomogenous multiplicity.
1. Resource usage as a non-functional error, meaning there are too many memory allocations or a mistake in choosing template parameters, leading to slow code.
1. Misuse of Julia features, like custom indexing, that leads to calling functions incorrectly.


## Addressing risk outside of testing

**Use a single notation.** I'll record notation in the documentation, as a separate document, and then translate each algorithm, from different sources, to use that notation. It's painstaking but important.

**Include testing code in the main library.** This is a way to make it easier to write tests when there is a lot of support code. By "testing code," I mean alternative versions of functions with which to compare other functions, so these aren't `@test` invocations bug supporting functions. I wish there were, in general, a way to designate part of a library as not autoloaded in a package, but we'll make do.


## Broad Test methods

* Comparison
  - Multiple implementations
    - Two functions with same input and output.
    - Implementation in another Julia library, like Dierckx.
    - Implementation in R or Python.
  - Multiple stages of derivation of the equations
    - Symbolic calculation
    - Polynomial-based calculation
    - B-spline-based calculation
  - Published values
    - Wikipedia
    - Gradshteyn and Ryzhik
* Theoretically tractable
  - Check B-spline invariants, such as the integral of the B-spline, its completeness over an interval.
  - Polynomial fit to a polynomially-generated set of data.
  - Construct B-spline, convert to polynomial, then convert back.
* Limiting forms
  - A non-uniform method on a uniform grid should agree with the uniform solution.
  - Increasing number of knots should approximate solution at a known rate.
* Scaling tests
  - Shift and scale of whole axis should shift and scale the functions on that axis.
* Numerical checks
  - Test with Float, BigFloat, Rational, Int.
  - Check for type stability.
  - Check that this can be auto-differentiated.
  - Look for divide-by-zero. That can happen here for close knots.
  - Matrix solution stability. Can use numerical checks. We'll work with recommended solutions from literature, but the numerical checks of stability are reassurance about which code is ready to use and which isn't.


## Test Harness Specifics

### Knot multiplicity translation

There are two ways to specify an axis. One repeats knots.

$(x_1, x_1, x_2, x_2\ldots,x_k,x_{k+1})$

The other uses an integer vector of multiplicities to indicate consecutive, identical knot values.

$(2,2,\ldots,1,1)$

I'll need to translate from one to the other and back.

### Polynomial versus Spline

The library rests on a few data structures, among which tests will need to translate.

1. The B-spline representation, which is a set of B-spline coefficients.
   a. An axis, `Float64[]` (or other type, but it's a vector. That's the point.).
   b. Coefficients, `Float64[]`.

2. The polynomial pieces representation, which is a set of polynomials defined on consecutive intervals.
   a. An axis, `Float64[]`.
   b. Polynomial constants, `Array{Float64}(coefficient order, interval count)`.

3. A polynomial representation of the B-splines. This isn't standard. It would be an array of polynomial pieces on each interval. If the order is `m`, then there are `m` B-splines at each interval.
   a. An axis, `Float64[]`.
   b. Polynomials for all B-splines on an interval, `Array{Float64}(coefficient order, B-spline index, interval count)`. This would be three-dimensional with size `(m, m, k+1)`, where there are `k+2` points on the axis defining `k+1` intervals.

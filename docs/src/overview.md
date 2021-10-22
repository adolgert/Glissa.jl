# Overview

This is a library for univariate splines on non-uniform axes. A spline defines a polynomial at each interval of an axis.

* Cubic splines, using Conte and de Boor's Elementary Numerical Analysis.

* B-spline bases, using Schumaker's _Spline Functions: Basic Theory,_ 3rd Ed. Schumaker also has a recent book on SplinePak for Matlab. It's not a book as much as documentation, but it's a good overview of how implementations fit together to do a job.

* Dierckx's book on spline functions.

I'm writing this library as a source of functions for contributing to the [Interpolations.jl](https://github.com/JuliaMath/Interpolations.jl) package, because that package currently supports uniform splines and doesn't have some traditional functionality for B-spline work. As a result, this library will focus on functionality of spline functions, which means lots of testing, and we do that by comparison. So here's what's inside this library.

1. Given B-spline coefficients, evaluate a polynomial on an axis. It doesn't build representations of the splines internally but calculates them on the fly. These algorithms come from Schumaker. `schumaker.jl`

2. Given B-spline coefficients, evaluate derivatives and integrals of polynomials on an axis. (Schumaker) `schumaker.jl`

3. Given B-spline coefficients, create a representation of the polynomial spline, as polynomial pieces. This means it generates a separate polynomial for each interval on the axis. Also from Schumaker. `schumaker.jl`

4. Interpolate data with cubic polynomials. This is an old-school polynomial function. This fit can also be done with splines, so we can use this for testing. (from Conte and de Boore's book). `cubic_spline.jl`

5. Construct B-spline values using recursive functions. These are old-school. Could use Direckx's implementation here. `divided_differences.jl`, `truncated.jl`.

6. Solve for the coefficients of a B-spline using a matrix of constraint equations. This isn't efficient for anything, but it's the definition of a B-spline. My derivation. `bsplinepolynomial.jl`.

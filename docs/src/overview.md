# Overview

This is a library for univariate splines on non-uniform axes. A spline defines a polynomial at each interval of an axis.

* Cubic splines, using Conte and de Boor's Elementary Numerical Analysis.

* B-spline bases, using Schumaker's _Spline Functions: Basic Theory,_ 3rd Ed.

You might think that this library will work a lot with polynomials, but the trick to splines is that you rarely represent the polynomials themselves. We'll represent them explicitly only for testing purposes. Instead, the code has the following parts:

1. Evaluate a spline at a point ``x``, given spline coefficients ``c_i``.

2. Evaluate the integral and derivatives of a spline given its coefficients.

3. Solve for spline coefficients, given values, derivatives, and boundary conditions.

There is a lot more that could be here, but that's where we'll start.

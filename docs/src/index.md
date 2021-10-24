```@meta
CurrentModule = Glissa
```

# Glissa

[Glissa](https://github.com/adolgert/Glissa.jl) is a Julia library that defines B-splines for non-uniform interpolation.

The [Interpolations.jl](https://github.com/JuliaMath/Interpolations.jl) package in Julia doesn't support non-uniform interpolation or interpolations with multiple knots (less smooth interpolation), so this package is a place to implement and test basic functions with which to expand functionality of Interpolations.jl.

There are many implementations of B-splines here.

1. **Iterative B-splines** - The most robust computation of a B-spline, implemented for non-uniform axes and repeated knots. This includes derivatives and integrals. [`bspline.jl`](https://github.com/adolgert/Glissa.jl/blob/main/src/bspline.jl).

2. **Recursive B-splines** - The most common way to define a B-spline is to use recursion. [`recursive.jl`](https://github.com/adolgert/Glissa.jl/blob/main/src/recursive.jl)

3. **Divided differences B-splines** - This is a function-on-a-function that is one way to define a B-spline. It's a theoretically-simple definition not normally used for computation, but Julia makes the computation nice enough. [`divided_differences.jl`](https://github.com/adolgert/Glissa.jl/blob/main/src/divided_differences.jl)

4. **Constraint-solved B-splines** - You can set up a set of equations on polynomial pieces, where these equations define the B-spline. This is an excellent check on other methods. [`bsplineconstriants.jl`](https://github.com/adolgert/Glissa.jl/blob/main/src/bsplineconstraints.jl)

There are also some other representations of splines, used for comparison.

* **Polynomial Pieces** - This represents a set of polynomials defined on contiguous intervals. Sometimes we convert splines into this space. [`piecewise.jl`](https://github.com/adolgert/Glissa.jl/blob/main/src/piecewise.jl)

* **Cubic spline fit** - Interpolate data with cubic polynomials. This is an old-school polynomial function. This fit can also be done with splines, so we can use this for testing. (from Conte and de Boore's book). [`cubic_spline.jl`](https://github.com/adolgert/Glissa.jl/blob/main/src/cubic_spline.jl)

* **Newton form** - A nice way to solve for a polynomial fit, from Conte and de Boor. [`newton_form.jl`](https://github.com/adolgert/Glissa.jl/blob/main/src/newton_form.jl)


The source material for this work is a set of standard implementations.

* Conte and de Boor's *Elementary Numerical Analysis.*

* Schumaker's _Spline Functions: Basic Theory,_ 3rd Ed. Schumaker also has a recent book on SplinePak for Matlab. It's not a book as much as documentation, but it's a good overview of how implementations fit together to do a job.

* Dierckx's book on spline functions.

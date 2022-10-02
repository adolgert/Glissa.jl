# Project Plan

## Goal

I'd like to expand the capabilities of [Interpolations.jl](https://github.com/JuliaMath/Interpolations.jl) so that it supports

a) Non-uniform interpolation

b) Varying smoothness of interpolant, also known as multiple knots.

## Scope

* This will have all of the B-spline functionality that a professional library should have, including derivatives, integrals, and conversion to polynomial pieces.
* This is not a standalone fitting library.
* This will not make its own user interface. It will be implementations of core algorithms.
* This will include statements of the math setup for each function but not derivations.


## Code Qualities

This code should be raw tools with which to build a different implementation, so that changes the usual code qualities.

1. Mutability will help these functions be adaptable to other uses, which is a main goal.
2. Testability is high on this list, because these routines are meant for use by another library.
3. Usability is a lower factor here than usual because this should work underneath another library's interface. The main usability concern is resource usage. So don't dictate the creation of buffers and do check for type stability.
4. Reliability isn't as much a problem because the problem domain is very controlled.
5. Security, this might make buffer overruns.


## Threats

* Redoing other work. - I did a survey of the spline libraries, and it looks like people made pieces of B-splines in order to serve particular purposes, but I don't see implementations that deal with the ugly non-uniform, multiple-knot stuff.

* Can't figure out how to use it in Interpolations.jl. - Ask for help here. And do a pre-read of the package to understand the structure as much as I can.

* Code is buggy. - I'll test the heck out of this.

* Code is incorrect. - Because there aren't a lot of implementations that include non-uniform, multiple-knot versions, there isn't a lot of comparison. Can ameliorate this by making multiple internal versions for testing. Also, use the well-known books: Schumaker, Dierckx, Conte and de Boor.

* Package loads slowly because it is larger than necessary. - I'm putting multiple versions into the same package, so that's excess code to compile. That's OK here because we'll pull out functionality, and that main functionality has very few dependencies.


## External Dependencies

We'll have lots of dependencies that are for testing.

* ForwardDiff.jl. This gives simple derivatives, so you don't have to think about it.
* Symbolics.jl - Such an excellent way to see whether two algorithms agree exactly. Pass a symbolic value in, and define the axis with rationals. Then the results should be exactly the same on the way out.
* RCall.jl, PyCall.jl - Compare against R and Python libraries.
* LinearAlgebra.jl - This is a core library. OK to use that liberally.
* Random, Distributions.jl - These generate sets of tests quickly.

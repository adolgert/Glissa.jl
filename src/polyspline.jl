"""
A general spline solver, the hard way, without basis functions.
Reading Chapter 4 on Polynomial Splines in Schumaker.
"""

"""
A polynomial spline defined without basis functions.
"""
mutable struct PolynomialSpline
    x::Vector  # The abcissa
    m::Vector  # Multiplicity of knots
    order::Int  # Order of polynomials, which is degree + 1. x^2 is order 3.
    c::Vector
end


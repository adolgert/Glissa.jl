module Glissa

# These are basic data structures.
include("axis_handling.jl")
include("piecewise.jl")
include("newton_form.jl")

# Definitions of B-splines.
include("truncated.jl")
include("divided_differences.jl")
include("bsplineconstraints.jl")
include("bspline.jl")

# Fitting with splines.
include("cubic_spline.jl")
include("schumaker.jl")

export PiecewisePolynomial
export cubic_spline
export evaluate!
export FreeSlope
export Monotonic

end

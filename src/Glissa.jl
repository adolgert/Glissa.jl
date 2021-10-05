module Glissa

include("divided_differences.jl")
include("truncated.jl")
include("newton_form.jl")
include("polyspline.jl")
include("bsplinepolynomial.jl")
include("bspline.jl")
include("cubic_spline.jl")
include("schumaker.jl")

export PolynomialSpline
export cubic_spline
export evaluate!
export FreeSlope
export Monotonic

end

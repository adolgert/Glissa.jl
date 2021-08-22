module Glissa

include("divided_differences.jl")
include("truncated.jl")
include("newton_form.jl")
include("bspline.jl")
include("cubic_spline.jl")

export CubicSpline
export cubic_spline_flat_endpoints
export cubic_spline
export evaluate!

end

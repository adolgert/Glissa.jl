using Glissa
using Test

@testset "Glissa.jl" begin
    include("test_axis_handling.jl")
    include("test_cubic_spline.jl")
    include("test_divided_differences.jl")
    include("test_recursive.jl")
    include("test_bsplineconstraints.jl")
    include("test_bspline.jl")
end

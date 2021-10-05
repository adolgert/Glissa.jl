using Glissa
using Test

@testset "Glissa.jl" begin
    include("test_cubic_spline.jl")
    include("test_bsplinepolynomial.jl")
    include("test_schumaker.jl")
end

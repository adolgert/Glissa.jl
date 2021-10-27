using Glissa
using Test

@testset "Glissa.jl" begin
    test_files = [
        "test_axis_handling.jl",
        "test_cubic_spline.jl",
        "test_divided_differences.jl",
        "test_recursive.jl",
        "test_bsplineconstraints.jl",
        "test_bspline.jl",
        "test_compare.jl"
    ]
    exclude_files = Vector{String}()
    missing_files = Vector{String}()
    exist_files = readdir(".")
    for filename in exist_files
        if endswith(filename, ".jl") && startswith(filename, "test_") &&
            filename ∉ test_files && filename ∉ exclude_files

            push!(missing_files, filename)
        end
    end
    if length(missing_files) > 0
        println("You forgot to include $(missing_files)")
        @test length(missing_files) == 0
    end
    for tf in test_files
        include(tf)
    end
end

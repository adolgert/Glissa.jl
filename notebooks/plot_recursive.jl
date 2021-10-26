using Glissa
using Plots
using Random

function plot_splines(uniques, multiples, order)
    bspline_cnt = sum(multiples) - order
    axis = Glissa.axis_multiples_to_repeats(uniques, multiples)
    dx = (axis[end] - axis[1]) / 100
    x = axis[1]:dx:axis[end]
    bspline_index = 1
    y = similar(x)
    plotout = nothing
    for bidx in 1:bspline_cnt
        for yi in 1:length(y)
            y[yi] = Glissa.bsplineq_recursive(axis, order, bidx, x[yi])
        end
        if bidx == 1
            plotout = plot(x, y, title = "Bspline recursive")
        else
            plotout = plot!(x, y)
        end
    end
    display(plotout)
end

plot_splines(1:10, ones(Int, 10), 3)
plot_splines([0, .4, 1], [2, 2, 2], 3)

rng = Random.MersenneTwister(432432)
axis, uniques, multiples, order = generate_random_polyspline_specification(rng, Float64)

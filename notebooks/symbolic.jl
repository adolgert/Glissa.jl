using Symbolics
using Revise
using Glissa


# This notebook finds B-spline polynomial coefficients by creating symbolic variables
# and asking a B-spline evaluator to compute from the symbolic variables.
#
# What are we going to match here? Take Wikipedia's order 3 splines:
# B_1 = x^2/2   0 ≤ x < 1
# B_2 = (-2x^2 + 6x - 3)/2  1 ≤ x < 2
# B_3 = (3-x)^2 / 2
#
# The wikipedia actually makes this hard b/c it writes B-splines using the same
# x=0 as the offset, but we compute with x-x_i where x_i is the knot to the left of x.
# So we will pick out the relevant pieces. We could, instead, compute this once and shift x.
@variables x

@assert typeof(x) == Num
order = 3  # polynomial order = degree + 1.
# This will hold the output equations. There will be `order` equations, but it needs
# an extra +1 at the end because it's used as a buffer inside the function.

# Argument is the x of the left-hand side of the interval.
function bspline_l(lefthandside::Int)
    y = UnitRange{Rational{Int}}(-order, 2*order)
    l = searchsorted(y, lefthandside).start
    N = zeros(Num, order + 1)
    Glissa.generate_normalized_bsplines!(N, y, l, order, x)
    for si in 1:order
        N[si] = simplify(N[si], expand = true)
    end
    N
end

# This is all B-splines that overlap the 0-to-1 interval.
# Really, this is all you would need to calculate.
N = bspline_l(0)
N[1:order]
# The last one is the left-hand interval of our B-spline.
B1 = N[3]

# Pick the index to get values from 1 to 2.
N = bspline_l(1)
B2 = N[2]

# Pick the index to get values from 2 to 3
N = bspline_l(2)
B3 = N[1]

[B1, B2, B3]
# julia> [B1, B2, B3]
# 3-element Vector{Num}:
#                      (1//2)*(x^2)
#         (3//1)*x - (3//2) - (x^2)
#  (9//2) + (1//2)*(x^2) - (3//1)*x

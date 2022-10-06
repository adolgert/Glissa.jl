using Logging


# What's a derivative? of sum(c_i x^(i-1), {i, 1, n})
# c[1] x^0 + c[2] x^1 + c[3] x^2 + c[4] x^3   # sum(c_i x^(i-1), {i, 1, n}), n=4.
#            1*c[2] x^0 2*c[3] x^1 + 3*c[4] x^2  # sum((i-1) c_i x^(i-2), {i, 2, n})
#                       2 c[3] x^0 + 2*3 c[4] x^1 # sum((i-2)(i-1) c_i x^(i-3), {i, 3, n})
# jth derivative: for j>0.
#   sum((factorial(i-1)/factorial(i-1-j) x^(i-1-j)), {i,1-j,n})
#   sum(prod(k, {k,i-j, i-1}) x^(i-1-j)), {i,1-j,n})

# A product from a to b.
@doc raw"""``a(a+1)(a+2)\cdots(b-1)b``"""
function dprod(a, b)
    x::typeof(a * b) = 1
    for i::typeof(a * b) = a:b
        x *= i
    end
    return x
end

# Places values in the matrix, at the given row, that mean the value of this polynomial
# is part of a constraint equation. Scale that value by the scale.
function valueat!(mat::AbstractArray{T}, axis, interval, row, side, scale, order) where {T}
  # This if-then is because equations are in Horner form (x-x_0) wher x_0 is the left side
  # of the interval in which this polynomial is defined.
    x = if side == :Left
        zero(T)
    else
        axis[interval + 1] - axis[interval]
    end
    for i in 1:order
        col = i + (interval - 1) * order
        mat[row, col] += x^(i - 1) * scale
    end
end

# j is derivative, so first derivative is j=1.
# interval is index of the span between axis points.
# row is the row of the matrix to add values to.
# side is whether this is the left or right side of the interval.
# scale multiplies values, so we can add or subtract, to form equations.
function derivat!(mat::AbstractArray{T}, axis, j, interval, row, side, scale, order) where {T}
    if j == 0
        valueat!(mat, axis, interval, row, side, scale, order)
        return nothing
    end
    x = if side == :Left
        zero(T)
    else
        axis[interval + 1] - axis[interval]
    end
    for i in (1 + j):order  # This drops terms b/c of the derivative.
        col = i + (interval - 1) * order
        mat[row, col] += dprod(i - j, i - 1) * x^(i - j - 1) * scale
    end
    return nothing
end

# Q(x) is normalized so that its integral is 1 / order
function integralat!(mat::AbstractArray{T}, axis, interval, row, order) where {T}
    @debug "$(length(axis)) $interval"
    for i in 1:order
        col = i + (interval - 1) * order
        mat[row, col] += (axis[interval + 1] - axis[interval])^i / i
    end
    return nothing
end


# N(x) is normalized so that the sum over the internal axis points is 1.
function pointsum!(mat::AbstractArray{T}, axis, interval, row, order) where {T}
    col = 1 + (interval - 1) * order
    mat[row, col] += 1
    return nothing
end

@doc raw"""
    bspline_by_matrix!(
      axis::AbstractArray{T},
      coeffs::AbstractArray{T},
      multiplicity::AbstractVector{Int},
      order, normalized
    ) where {T}

Calculate coefficients of a b-spline of order m, degree m-1
by sending an axis of length m+1 and a coefficient matrix of size mxm,
into which are written columns of coefficients for the polynomial of each interval.
If you want a b-spline on a subset of an axis, use an array view for inputs.
This is the least efficient way to calculate these coefficients, but it makes
explicit the continuity conditions that lead to a solution.
`i` is the index of the leftmost point, counting points by their multiplicity,
so a point with multiplicity of 3 is three points.
"""
function bspline_by_matrix!(
  axis::AbstractArray{T}, coeffs::AbstractArray{T}, multiplicity::AbstractVector{Int},
  order, normalized
  ) where {T}

    degree = order - 1
    multiples = sum(multiplicity .- 1)
    interval_cnt = order - multiples  # The number of intervals.
    @assert interval_cnt == length(axis) - 1
    @assert order == length(axis) - 1 + multiples
    @assert length(axis) == length(multiplicity)
    @assert order == size(coeffs, 1)
    @assert interval_cnt == size(coeffs, 2)
    mat = zeros(T, interval_cnt * order, interval_cnt * order)
    coeffs .= zero(T)

    # At the left-hand side
    row_idx = 1
    for deridx in 0:(degree - multiplicity[1])
        # Left side of first interval.
        derivat!(mat, axis, deridx, 1, row_idx, :Left, 1, order)
        row_idx += 1
    end

    # Internal nodes
    for k in 1:(interval_cnt - 1)
        # At each join, the derivatives match.
        # k + 1 because it's multiplicity of the vertex, which is 1 + interval index.
        for deridx in 0:(degree - multiplicity[k + 1])
            # Here, k is the index of the interval.
            derivat!(mat, axis, deridx, k, row_idx, :Right, -1, order)
            derivat!(mat, axis, deridx, k + 1, row_idx, :Left, 1, order)
            row_idx += 1
        end
    end

    # At the right-hand side
    for deridx in 0:(degree - multiplicity[interval_cnt + 1])
        # Right side of last interval.
        derivat!(mat, axis, deridx, interval_cnt, row_idx, :Right, 1, order)
        row_idx += 1
    end

    # The last constraint is that the integral over the b-spline is one.
    for int_idx in 1:interval_cnt
        integralat!(mat, axis, int_idx, row_idx, order)
    end
    rhs = zeros(T, interval_cnt * order)
    if normalized == :Q
        rhs[end] = one(T) / T(order)  # The integral is 1/m. for Q(x).
    elseif normalized == :N
        # Normalization is (y_m - y_1) / m.
        rhs[end] = (axis[end] - axis[1]) / T(order)
    else
      error("Normalization not recognized")
    end
    # The matrix must be square in order for it to have an inverse.
    @assert row_idx == size(mat, 1)
    # The matrix will be nonsingular iff it has a positive entry at every diagonal.
    # print("bspline_by_matrix! construction ")
    # display("text/plain", mat)
    # println()
    x = mat \ rhs
    for interval_idx in 1:interval_cnt
        for kidx in 1:order
            coeffs[kidx, interval_idx] = x[(interval_idx - 1) * order + kidx]
        end
    end
end


@doc raw"""
    reduce_axis(multiplicity::AbstractVector, order, i)

This extracts from a multiplicity vector the support for the `i`th B-spline.
There is an axis where `length(axis)=length(multiplicity)`. This axis will, in general,
support multiple B-splines. This function returns a sub-range of that axis
and a new multiplicity vector, the same length as that sub-range. It will be the exact
support for the `i`th B-spline of order `order`.
"""
function reduce_axis(multiplicity::AbstractVector, order, i)
    cumulative = cumsum(multiplicity)
    left_index = searchsortedfirst(cumulative, i)
    right_index = searchsortedfirst(cumulative, i + order)
    multiplicityprime = copy(multiplicity[left_index:right_index])

    multiplicityprime[1] = cumulative[left_index] - i + 1
    multiplicityprime[end] -= sum(multiplicityprime) - (order + 1)

    ((left_index, right_index), multiplicityprime)
end



@doc raw"""
    polyspline_constraints!(axis::AbstractArray{T}, multiplicity::AbstractVector{Int}, order) where {T}

This reads a list of multiplicities and converts that into a matrix that represents
the constraints indicated by those multiplicites. The resulting matrix defines the space
of the given polynomial spline on this axis.
"""
function polyspline_constraints!(axis::AbstractArray{T}, multiplicity::AbstractVector{Int}, order) where {T}

    degree = order - 1
    @assert length(axis) == length(multiplicity)
    interval_cnt = length(axis) - 1
    boundary_conditions_cnt = 2 * order - multiplicity[1] - multiplicity[end]
    knot_conditions_cnt = sum(order .- multiplicity[2:(length(multiplicity) - 1)])
    mat = zeros(T, boundary_conditions_cnt + knot_conditions_cnt, interval_cnt * order)

    # At the left-hand side
    row_idx = 1
    for deridx in 0:(degree - multiplicity[1])
        # Left side of first interval.
        derivat!(mat, axis, deridx, 1, row_idx, :Left, 1, order)
        row_idx += 1
    end

    # Internal nodes
    for k in 1:(interval_cnt - 1)
        # At each join, the derivatives match.
        # k + 1 because it's multiplicity of the vertex, which is 1 + interval index.
        for deridx in 0:(degree - multiplicity[k + 1])
            # Here, k is the index of the interval.
            derivat!(mat, axis, deridx, k, row_idx, :Right, -1, order)
            derivat!(mat, axis, deridx, k + 1, row_idx, :Left, 1, order)
            row_idx += 1
        end
    end

    # At the right-hand side
    for deridx in 0:(degree - multiplicity[interval_cnt + 1])
        # Right side of last interval.
        derivat!(mat, axis, deridx, interval_cnt, row_idx, :Right, 1, order)
        row_idx += 1
    end

    @assert row_idx == size(mat, 1) + 1
    mat
end

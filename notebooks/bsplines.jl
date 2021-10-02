# This notebook solves for B-splines of order 2 by constructing a matrix
# of constraints and solving that matrix for coefficients of each part of the
# polynomial.

# What's a derivative? of sum(c_i x^(i-1), {i, 1, n})
# c[1] x^0 + c[2] x^1 + c[3] x^2 + c[4] x^3   # sum(c_i x^(i-1), {i, 1, n}), n=4.
#            1*c[2] x^0 2*c[3] x^1 + 3*c[4] x^2  # sum((i-1) c_i x^(i-2), {i, 2, n})
#                       2 c[3] x^0 + 2*3 c[4] x^1 # sum((i-2)(i-1) c_i x^(i-3), {i, 3, n})
# jth derivative: for j>0.
#   sum((factorial(i-1)/factorial(i-1-j) x^(i-1-j)), {i,1-j,n})
#   sum(prod(k, {k,i-j, i-1}) x^(i-1-j)), {i,1-j,n})

# A product from a to b.
function dprod(a, b)
  x::typeof(a*b) = 1
  for i::typeof(a*b) = a:b
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
    mat[row, col] = x^(i - 1) * scale
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
    mat[row, col] = dprod(i-j, i-1) * x^(i - j - 1) * scale
  end
  return nothing
end

# Q(x) is normalized so that its integral is 1 / order
function integralat!(mat::AbstractArray{T}, axis, interval, row, order) where {T}
  order = length(axis) - 1
  x = axis[interval + 1] - axis[interval]
  for i in 1:order
    col = i + (interval - 1) * order
    mat[row, col] = (axis[interval + 1]^i - axis[interval]^i) / i
  end
  return nothing
end


# N(x) is normalized so that the sum over the internal axis points is 1.
function pointsum!(mat::AbstractArray{T}, axis, interval, row, order) where {T}
  x = axis[interval + 1] - axis[interval]
  col = 1 + (interval - 1) * order
  mat[row, col] = 1
  return nothing
end


# Calculate coefficients of a b-spline of order m, degree m-1
# by sending an axis of length m+1 and a coefficient matrix of size mxm,
# into which are written columns of coefficients for the polynomial of each interval.
# If you want a b-spline on a subset of an axis, use an array view.
# This is the least efficient way to calculate these coefficients, but it makes
# explicit the continuity conditions that lead to a solution.
function bspline_by_matrix!(
  axis::AbstractArray{T}, coeffs::AbstractArray{T}, multiplicity::AbstractVector{Int},
  order, normalized
  ) where {T}

  degree = order - 1
  multiples = sum(multiplicity .- 1)
  vcnt = order - multiples
  @assert order == length(axis) - 1 + multiples
  @assert vcnt == size(coeffs, 2)
  @assert vcnt == length(multiplicity)
  mat = zeros(T, vcnt * order, vcnt * order)

  # At the left-hand side
  row_idx = 1
  for deridx in 0:(degree - multiplicity[1])
    derivat!(mat, axis, deridx, 1, row_idx, :Left, 1, order)
    row_idx += 1
  end

  # Internal nodes
  for k in 1:(degree - multiples)
    # At each join, the derivatives match.
    for deridx in 0:(degree - multiplicity[k])
      derivat!(mat, axis, deridx, k, row_idx, :Right, -1, order)
      derivat!(mat, axis, deridx, k + 1, row_idx, :Left, 1, order)
      row_idx += 1
    end
  end

  # At the right-hand side
  for deridx in 0:(degree - multiplicity[vcnt])
    derivat!(mat, axis, deridx, vcnt, row_idx, :Right, 1, order)
    row_idx += 1
  end

  if normalized == :Q
    # The last constraint is that the integral over the b-spline is one.
    for int_idx in 1:order
      integralat!(mat, axis, int_idx, row_idx, order)
    end
    rhs = zeros(T, vcnt * order)
    rhs[end] = one(T) / T(order)  # The integral is 1/m. for Q(x).
  elseif normalized == :N
    for int_idx in 2:(order - 1)
      pointsum!(mat, axis, int_idx, row_idx, order)
    end
    rhs = zeros(T, vcnt * order)
    rhs[end] = one(T)   # The sum is 1 for N(x).
  else
    error("Normalization not recognized")
  end
  @assert row_idx == size(mat, 1)
  x = mat \ rhs
  for interval_idx in 1:vcnt
    for kidx in 1:order
      coeffs[kidx, interval_idx] = x[(interval_idx - 1) * order + kidx]
    end
  end
end


# Use rationals so that we can compare with theoretical calculations easily.
T = Rational{Int64}
deg = 2  # degree of polynomial.
ord = deg + 1  # order of polynomial.
axis = collect(zero(T):T(ord))
coeffs = zeros(T, ord, ord)
multiplicity = ones(Int, ord)
bspline_by_matrix!(axis, coeffs, multiplicity, ord, :N)

# We can compare that with Wikipedia's order 3 cardinal B-spline using symbolic math.
using Symbolics
@variables x
# First equation clearly matches.
sum(coeffs[1:3] .* [1, x, x^2])
# Second equation, if we apply the horner shift and expand it, matches.
simplify(sum(coeffs[4:6] .* [1, (x-1), (x-1)^2]), expand=true)
simplify((-2x^2+6x-3)/2, expand = true)
# So does third equation.
simplify(sum(coeffs[7:9] .* [1, (x-2), (x-2)^2]), expand=true)
simplify((3-x)^2/2, expand = true)

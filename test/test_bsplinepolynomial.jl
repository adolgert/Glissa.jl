using Glissa
using Random
using Distributions

# test axis reduction with reduce_axis
rng = Random.MersenneTwister(2094372)
max_order = 5
for i in 1:10000
  order = rand(rng, 1:max_order)
  vertex_cnt = rand(rng, (order + 1):(2*order + 2))
  multiplicity = Int[]
  vertex_idx = vertex_cnt
  while vertex_idx > 0
    if order > 2 && rand(rng, 1:5) == 1
      mult = min(rand(rng, 1:(order - 1)), vertex_idx)
      append!(multiplicity, mult)
      vertex_idx -= mult
    else
      append!(multiplicity, 1)
      vertex_idx -= 1
    end
  end
  @assert sum(multiplicity) == vertex_cnt
  println("ord $order mult $multiplicity")
  for i in 1:(vertex_cnt - order)
    ((l, r), m) = Glissa.reduce_axis(multiplicity, order, i)
    println("i $i l $l r $r m $m")
    @assert sum(m) == order + 1
    @assert r - l + 1 == length(m)
    for j in l:r
      @assert m[j - l + 1] ≤ multiplicity[j]
    end
  end
end


struct PolyBSpline{X,T} <: PolynomialSpline{X,T}
  τ::AbstractVector{X}  # The abcissa
  c::AbstractMatrix{T}
  bounds::UnitRange{Int}
end

function evaluate(totalaxis::AbstractVector, bsplines::AbstractVector, x)
  i = searchsortedlast(totalaxis, x)
  total::eltype(bsplines[1]) = 0
  for j in 1:length(bsplines)
    if i ≥ bsplines[j].bounds.start && i < bsplines[j].bounds.stop
      total += bsplines[j](x)
    end
  end
  total
end



# This tests solution of the constituent equations along an axis.
rng = Random.MersenneTwister(19472234)
max_order = 5
for i in 1:5
  order = rand(rng, 1:max_order)
  vertex_cnt = rand(rng, (order + 1):(2*order + 2))
  multiplicity = Int[]
  vertex_idx = vertex_cnt
  while vertex_idx > 0
    if order > 2 && rand(rng, 1:5) == 1
      mult = min(rand(rng, 1:(order - 1)), vertex_idx)
      append!(multiplicity, mult)
      vertex_idx -= mult
    else
      append!(multiplicity, 1)
      vertex_idx -= 1
    end
  end
  @assert sum(multiplicity) == vertex_cnt
  println("ord $order mult $multiplicity")
  T = Float64
  a::Float64 = 0.0
  b::Float64 = length(multiplicity) - 1
  axis = collect(UnitRange{T}(a, b))
  splines = Vector{PolyBSpline{T,T}}(undef, vertex_cnt - order)
  for i in 1:(vertex_cnt - order)
    ((l, r), m) = Glissa.reduce_axis(multiplicity, order, i)
    coeffs = zeros(T, order, r - l + 1)
    println("$(length(axis)) $(size(coeffs)) $m $order")
    Glissa.bspline_by_matrix!(axis[l:r], coeffs, m, order, :Q)
    # Glissa.bspline_by_matrix!(axis[l:r], coeffs, m, order, :N)
    splines[i] = PolyBSpline{T,T}(axis[l:r], coeffs, UnitRange{Int}(l, r))
  end

  for check_idx in 1:100
    x = rand(rng, Uniform(a, b))
    y = evaluate(axis, splines, x)
    @test y ≈ 1.0
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

T = Rational{Int64}
deg = 2  # degree of polynomial.
ord = deg + 1  # order of polynomial.
axis = collect(zero(T):(T(ord) - 1))
coeffs = zeros(T, ord, ord - 1)
multiplicity = ones(Int, ord - 1)
multiplicity[2] = 2
bspline_by_matrix!(axis, coeffs, multiplicity, ord, :N)

coeffs
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

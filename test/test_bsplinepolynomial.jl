using Glissa
using Random
using Distributions
using Symbolics
using Test

@testset "Can extract subset of an axis for bsplines" begin
  # test axis reduction with reduce_axis
  rng = Random.MersenneTwister(2094372)
  max_order = 5
  for trial_idx in 1:10
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
    @test sum(multiplicity) == vertex_cnt
    println("ord $order mult $multiplicity")
    for i in 1:(vertex_cnt - order)
      ((l, r), m) = Glissa.reduce_axis(multiplicity, order, i)
      println("i $i l $l r $r m $m")
      @test sum(m) == order + 1
      @test r - l + 1 == length(m)
      for j in l:r
        @test m[j - l + 1] ≤ multiplicity[j]
      end
    end
  end
end

struct PolyBSpline{X,T} <: PiecewisePolynomial{X,T}
  τ::AbstractVector{X}  # The abcissa
  c::AbstractMatrix{T}
  bounds::UnitRange{Int}
end


function evaluate(totalaxis::AbstractVector, bsplines::AbstractVector, x)
  i = searchsortedlast(totalaxis, x)
  total::eltype(bsplines[1]) = 0
  for j in 1:length(bsplines)
    # i is the index of the polynomial piece. The bounds are axis vertices, so 1 + pieces.
    if i ≥ bsplines[j].bounds.start && i < bsplines[j].bounds.stop
      total += bsplines[j](x)
    end
  end
  total
end


@testset "B-splines are greater than zero on their domain" begin
  # This tests solution of the constituent equations along an axis.
  rng = Random.MersenneTwister(19472234)
  max_order = 5
  for trial_idx in 1:5
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
    @test sum(multiplicity) == vertex_cnt
    println("ord $order mult $multiplicity")
    T = Float64
    a::Float64 = 0.0
    b::Float64 = length(multiplicity) - 1
    axis = collect(UnitRange{T}(a, b))
    splines = Vector{PolyBSpline{T,T}}(undef, vertex_cnt - order)
    for i in 1:(vertex_cnt - order)
      ((l, r), m) = Glissa.reduce_axis(multiplicity, order, i)
      coeffs = zeros(T, order, r - l) # A set of coefficients for each interval on the axis.
      println("$(length(axis)) $(size(coeffs)) $m $order")
      # Using the Q splines, which are normed to 1/order.
      Glissa.bspline_by_matrix!(axis[l:r], coeffs, m, order, :Q)
      splines[i] = PolyBSpline{T,T}(axis[l:r], coeffs, UnitRange{Int}(l, r))
      for check_x in (axis[l] + 1e-3):0.1:(axis[r] - 1e-3)
        @test splines[i](check_x) > 0
      end
    end
  end
end


@testset "Normal B-splines sum to one" begin
  # This tests solution of the constituent equations along an axis.
  rng = Random.MersenneTwister(19472234)
  max_order = 5
  for trial_idx in 1:5
    order = rand(rng, 1:max_order)
    piece_cnt = rand(rng, 1:(2*order + 2))
    multiplicity = zeros(Int, piece_cnt + 1)
    # It's here that we pad the edges with full multiplicity so that
    # the whole domain will sum to one over all N_i(x).
    multiplicity[1] = order
    multiplicity[piece_cnt + 1] = order
    for midx in 2:piece_cnt
      if order > 2 && rand(rng, 1:5) == 1
        mult = rand(rng, 1:order)
        multiplicity[midx] = mult
      else
        multiplicity[midx] = 1
      end
    end
    println("ord $order mult $multiplicity")
    T = Float64
    a::Float64 = 0.0
    b::Float64 = length(multiplicity) - 1
    axis = collect(UnitRange{T}(a, b))
    spline_cnt = sum(multiplicity) - order
    splines = Vector{PolyBSpline{T,T}}(undef, spline_cnt)
    for i in 1:spline_cnt
      ((l, r), m) = Glissa.reduce_axis(multiplicity, order, i)
      coeffs = zeros(T, order, r - l)
      println("$(length(axis)) $(size(coeffs)) $m $order")
      Glissa.bspline_by_matrix!(axis[l:r], coeffs, m, order, :N)
      splines[i] = PolyBSpline{T,T}(axis[l:r], coeffs, UnitRange{Int}(l, r))
    end

    for check_idx in 1:5
      x = rand(rng, Uniform(a, b - 1e-6))
      y = evaluate(axis, splines, x)
      @test y ≈ 1.0
    end
  end
end


@testset "B-splines match the known ones" begin

  T = Rational{Int64}
  deg = 2  # degree of polynomial.
  ord = deg + 1  # order of polynomial.
  axis = collect(zero(T):T(ord))
  coeffs = zeros(T, ord, ord)
  multiplicity = ones(Int, ord + 1)
  Glissa.bspline_by_matrix!(axis, coeffs, multiplicity, ord, :N)

  # We can compare that with Wikipedia's order 3 cardinal B-spline using symbolic math.
  Symbolics.@variables x
  # First equation clearly matches.
  @test simplify(sum(coeffs[1:3] .* [1, x, x^2]) == x^2 / 2) == true
  # Second equation, if we apply the horner shift and expand it, matches.
  c = simplify(sum(coeffs[4:6] .* [1, (x-1), (x-1)^2]), expand=true)
  d = simplify((-2x^2+6x-3)/2, expand = true)
  @test simplify(c == d) == true
  # So does third equation.
  a = simplify(sum(coeffs[7:9] .* [1, (x-2), (x-2)^2]), expand=true)
  b = simplify((3-x)^2/2, expand = true)
  @test simplify(a == b) == true
  
end

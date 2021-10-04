using Glissa
using Random
using Distributions
using Plots

struct PolyBSpline{X,T} <: PolynomialSpline{X,T}
  τ::AbstractVector{X}  # The abcissa
  c::AbstractMatrix{T}
  bounds::UnitRange{Int}
end

function evaluate(totalaxis::AbstractVector, bsplines::AbstractVector, x)
  i = searchsortedlast(totalaxis, x)
  total::eltype(bsplines[1]) = 0
  for j in 1:length(bsplines)
    if i in bsplines[j].bounds
      total += bsplines[j](x)
    end
  end
  total
end

function plotsplines(splines, i, append = false)
	xp = splines[i].τ[1]:0.1:splines[i].τ[end]
	yp = splines[i].(xp)
	if append
		println("plot $i append")
		display(plot!(xp, yp))
	else
		println("plot $i new")
		display(plot(xp, yp))
	end
end

function plotallsplines(splines)
	n = length(splines)
	plotsplines(splines, 1)
	for i in 2:n
		plotsplines(splines, i, true)
	end
end

order = 3
T = Float64
axis = T[0, 1, 2, 3, 4]
multiplicity = [order, 2, 2, 2, order]
@assert length(multiplicity) == length(axis)
vcnt = sum(multiplicity)
a = 0
b = 5

splines = Vector{PolyBSpline{T,T}}(undef, sum(multiplicity) - order)
for i in 1:(sum(multiplicity) - order)
	((l, r), m) = Glissa.reduce_axis(multiplicity, order, i)
	coeffs = zeros(T, order, r - l + 1)
	println("$m")
	Glissa.bspline_by_matrix!(axis[l:r], coeffs, m, order, :Q)
	splines[i] = PolyBSpline{T,T}(axis[l:r], coeffs, UnitRange{Int}(l, r))
end


length(splines)
plotsplines(splines, 3, true)
plotallsplines(splines)

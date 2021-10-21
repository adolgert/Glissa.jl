@doc raw"""
The Newton form for an interpolating polynomial represents ``f(x)`` at
points ``\tau``` by a sum.

``\sum_{i=0}^n f[x_0,\ldots,x_i]\prod_{j=0}^{i-1}(x-x_j)

Here, the `a` are the divided differences and ``\tau`` are the `x[j]`.
"""
struct NewtonForm{T}
    τ::Vector
    a::Vector{T}
end


"""
Construct a NewtonForm. From Conte and deBoor 1980, Eqn. 4.8.
"""
function NewtonForm(τ::Vector, f::Vector{T}) where {T <: Real}
    k = len(τ)
    a = zeros(T)
    dd = DividedDifference{T}()
    for j = 1:k
        a[j] = divided_difference(1, j, τ, f, dd)
    end
    NewtonForm{T}(τ, a)
end


function (nf::NewtonForm{T})(x) where {T}
    k = len(nf.a)
    total = zero(T)
    for j = k:-1:2
        total += a[j]
        total *= (x - τ[j - 1])
    end
    total += a[1]
    total
end


@doc raw"""
Given an axis `τ` and values `f` on that axis, return a polynomial `y` such that
``y(τ[i]) = f[i]``.
Construct a polynomial that interpolates a function, using polynomials with zeroes.
This should match the constructor above. It comes from deBoor equation 4.2.
"""
function single_polynomial(τ::Vector{X}, f::Vector{T}) where {T <: Real, X<: Real}
    a = zeros(T, length(f))
    for i = 1:length(f)
        total = one(T)
        for l1 = 1:(i - 1)
            total *= 1 / (τ[l1] - τ[i])
        end
        for l2 = (i + 2):length(f)
            total *= 1 / (τ[l2] - τ[i])
        end
        a[i] = f[i] * total
    end
    NewtonForm{T}(τ, a)
end

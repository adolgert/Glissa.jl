# This is recursive definitions of the B-spline.

@doc raw"""
    bsplineq(axis, order, index, x)

This is the ``Q_i^m(x)`` B-spline, which is normalized to ``1/m``. Here is is defined by
recursion.
"""
function bsplineq_recursive(y, m, i, x)
    if !(y[i] ≤ x < y[i+m])
        return zero(eltype(y))
    end

    if m == 1
        1/(y[i+1] - y[i])

    # repeated right knots: y[i] < y[i+1] == y[i+2] ... y[i+m]
    elseif length(y) > 2 && y[(i+2):(i+m)] .== y[i + 1]
        (x - y[i])^(m - 1) / (y[i+m] - y[i])^(m - 1)

    # repeated left knots: y[i] == y[i+1] == y[i+2] ... y[i+m-1] < y[i+m]
    elseif length(y) > 2 && y[(i+1):(i+m-1)] .== y[i]
        (y[i+m] - x)^(m - 1) / (y[i+m] - y[i])

    else # recurse
        ((x - y[i]) * bspline(y, m - 1, i, x) + (y[i + m] - x) * bspline(y, m - 1, i + 1, x)) /
            (y[i + m] - y[i])
    end
end


@doc raw"""
    bsplineq(axis, order, index, x, nth_derivative)

This is the ``Q_i^m(x)`` B-spline, which is normalized to ``1/m``. Here is is defined by
recursion.
"""
function bsplineq_recursive(y, m, i, x, n)
    if n == 0
        bsplineq_recursive(y, m, i, x)
    elseif n == 1
        (m - 1) * (bsplineq_recursive(y, m - 1, i, x) - bsplineq_recursive(y, m - 1, i + 1, x)) /
            (y[i + m] - y[i])
    else
        (m - 1) * (bsplineq_recursive(y, m - 1, i, x, n - 1) - bsplineq_recursive(y, m - 1, i + 1, x, n - 1)) /
            (y[i + m] - y[i])
    end
end

function bsplinen_recursive(y, m, i, x)
    (y[i + m] - y[i]) * bsplineq_recursive(y, m, i, x)
end

function bsplinen_recursive(y, m, i, x, n)
    (y[i + m] - y[i]) * bsplineq_recursive(y, m, i, x, n)
end

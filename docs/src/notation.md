# Notation

This defines terms in order to help read the code.

## Polynomial

A **polynomial** is a function of $x$, relative to some $x_j$, with constants $c_i$, of the form

$p_j(x) = $c_1 + c_2 (x - x_j) + c_3 (x - x_j)^2 = \sum_{i=1}^{i=d+1}c_i(x-x_j)^{i-1}$

The polynomial above has a **degree** of $d=2$ because $(x-x_1)^2$ is the largest power. It has **order** $m=3$ because there are 3 constants $(c_1, c_2, c_3)$. In Julia, you can evaluate a polynomial with [`evalpoly`](https://docs.julialang.org/en/v1/base/math/#Base.Math.evalpoly).

```julia
x = 5.5
x1 = 5.0
c = rand(3)
y = evalpoly(3, c)
```

## Piecewise Polynomial

Separate the $x$-axis into $(x_1, x_2, x_3, x_4... x_{k+2})$ so that there are $k + 1$ intervals, where **interval** $j$ is the axis between $x_j$ and $x_{j+1}$, not including the endpoint $x_{j+1}$. Now define a polynomial on each interval. Each of these is a **polynomial piece**. The $x_2\ldots x_{k+1}$-values, which excludes the endpoints, are called **knots.**

We need two indices on the constants, so they are now $c_{ij}$, where the row is the polynomial constant and the column is the polynomial piece. In addition, we assume all polynomial pieces have the same degree.

## Polynomial Spline

The piecewise polynomial above makes no guarantees that neighboring intervals will be continuous. If the value of the polynomial in one interval is continuous with the value of the neighboring polynomial, it's called $C^0$ smooth. If, in addition to being continuous, the first derivatives match, then it's $C^1$ smooth. For matching second derivatives, that's $C^2$ smooth.

A **polynomial spline with simple knots** is a piecewise polynomial, of degree $d$, which is $C^{d-1}$ smooth at each knot. It couldn't be any more smooth, because if we say, for instance, that two second-degree polynomials have the same values and derivatives, all the way up to second derivatives, then they are the same polynomial because they have to have the same polynomial constants.

A general **polynomial spline** guarantees that you don't have a situation where the derivatives agree across a knot, but the values aren't continuous. That may seem weird, but it's useful for some kinds of integration. Because of this restriction, we can think of each knot as having its own continuity condition. For instance, for a quadratic spline, each knot could be $C^0$ or $C^1$. For a quadratic spline, each knot could be $C^0$, $C^1$, $C^2$, or $C^3$.

The **multiplicity** of a knot indicates its smoothness. For a degree $d$ polynomial, a multiplicity $m_j=1$ know has maximal smoothness, which is $C^{d-1}$. As multiplicity increases, smoothness across the knot decreases. The maximum multiplicity is $m_j=m$, where $m=d+1$ is the order of the polynomial. When $j_j=m$, then the two neighboring polynomial pieces are disconnected. Note that endpoints aren't knots, just places where intervals meet. In code, an array `m[i]` is a multiplicity vector, but an integer `m` is the order.

There is another way to represent knot multiplicity. When listing the knots, repeat values for the knots of higher multiplicity, so $(x_1, x_1, x_1, x_2, x_2, x_3, x_3 x_4, x_4, x_4)$. Some algorithms rely on neighboring knots sometimes being equal.

There is some sense to using the word, multiplicity, to denote smoothness of splines. If three neighboring intervals are $C^2$ smooth, and you drag the knots together, shrinking the middle interval, the resulting neighboring intervals will only be $C^1$ smooth. Overlapping knots are equivalent to reduced continuity.

## Boundary conditions

At the left-most $x_1$ and the right-most $x_{k+1}$, we can define boundary conditions. A **boundary condition** is a specific value for the $j$-th derivative of the polynomial at that boundary. We will discuss so-called natural boundary conditions, meaning that we require the polynomial value to be 0 at the boundary, or its value and first derivative both to be 0, or its value and first two derivatives to be 0, and so on.

For a degree $d=2$ polynomial, we can require the value to be zero and the first derivative to be zero. That's up to 2 boundary conditions. If we required more, it would set the polynomial to be zero within the interval neighboring the boundary.

A boundary condition at $x_1$ or $x_{k+1}$ is kind of like a continuity condition with some polynomial piece that's outside the axis of the polynomial spline. Because of this similarity, we can represent boundary conditions by 1) defining values for a neighboring polynomial on each side and 2) assigning a multiplicity to the endpoints the same way we assign multiplicity to the knots. In this case, a multiplicity of 1 means that, at the endpoint, the polynomial piece is $C^{d-1}$ continuous with the polynomial outside. A multiplicity of $d+1$ means there is no continuity, so there are no boundary conditions at that endpoint.

| Degree | Boundary Multiplicity | Continuity with Zero       |
|--------|-----------------------|----------------------------|
| 2      | 3                     | not continuous             |
| 2      | 2                     | $C^0$ continuous           |
| 2      | 1                     | $C^1$, once differentiable |
| 3      | 4                     | not continuous             |
| 3      | 3                     | $C^0$                      |
| 3      | 2                     | $C^1$                      |
| 3      | 1                     | $C^2$ twice differentiable |

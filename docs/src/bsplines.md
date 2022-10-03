# B-splines

## Introduction

The [Wikipedia B-spline](https://en.wikipedia.org/wiki/B-spline) article tells you how to calculate them but doesn't describe what defines a B-spline. That definition will be useful for deriving, understanding, and testing algorithms, so let's cover it here.

### Compare with polynomials for interpolation

The [Notation](@ref) section described polynomial splines as a set of polynomials defined on neighboring intervals. For computational purposes, we don't usually represent polynomial splines as a list of polynomial constants. Instead, we use B-splines, because they are a more compact representation and offer advantages for computations using polynomial splines.

On an axis with points at `(0, 1, 2, 3, 4)`, with free boundary conditions, there are 6 degree-2 B-splines.

![Order 3 Splines](order3splines.png)

Each B-spline is a function that's defined over the whole interval, but, at the same time, it's only non-zero over at-most $d$ contiguous intervals.

Imagine that our goal is to fit data with a polynomial spline. We have points to fit, and they are between some $x_1$ and $x_{k+2}$, where we call the endpoint $k+2$ so that the count of internal knots is $k$. The first step is to divide the axis into intervals. At each knot where intervals meet, we get to decide, beforehand, whether the fit polynomial should be $C^0$-smooth (continuous), $C^1$-smooth (once differentiable), or as smooth as possible $C^{d-1}$-smooth. We can choose this knot-by-knot, if we want.

Now how would we do this fit with polynomials? For each interval, there is a separate polynomial $y_j$.

$y_j(x) = \sum_{i=1}^{i=m} c_{ij} (x-x_j)^{i-1}$

We could fit the points in order to find the $c_{ij}$ constants. For $k+1$ intervals on $k$ knots, there would be $(d+1)(k+1)$ constants to choose. When solving for the constants, we need to include equations that assert smoothness across the knots.

On the other hand, B-splines are functions that are designed, from their definition, to obey smoothness constraints. If constraints are built into the B-spline, we don't need to include equations to assert those same constraints. There are also fewer constants to choose becasue the final fit will be a sum of the B-splines.

$y(x) = \sum_j c_j B_j(x).$

That means there are exactly as many constants $c_j$ as we have degrees of freedom in the problem. If we are fitting with maximally-smooth polynomials, then there are $k+2d$ B-splines, total. And, as a neat trick, we don't need to store the polynomials that define the $B_j(x)$. They are quick to calculate from the knots.

**Prefiltering** is a name for finding the $c_j$ values that make a set of B-splines interpolate a set of points. You're taking data points $(x_p, y_p)$ and solving for the $c_j$ that make a B-spline-defined polynomial that matches them.


## B-spline definition from its constraints

B-splines are a basis set for polynomial splines on a particular axis, $x_1 < x_2 < x_3 \ldots < x_{k+2}$, with a given set of natural boundary conditions. Because they are a basis set, every allowable polynomial spline with those boundary conditions can be written as a sum of the B-splines for that axis. The set of B-splines reduces the amount of storage needed for a polynomial spline and reduces the computation needed to fit points with a polynomial spline. Now, how could we calculate the polynomial representation of any one of the B-splines for a given axis?

Let's compute the polynmial coefficients $c_{ij}$ that define a single B-spline on an axis. Start by picking a left-most point, some $x_j$. Then the B-spline will be zero to the left, positive for some number of intervals, say $p$ intervals. And then it will be zero to the right of that. We can count the unknown values as $p$ intervals of polynomials of degree $d$, so there are $p(d+1)$ unknown values. Oddly, we don't know the number of intervals for this B-spline yet, but continuity conditions will tell us what $p$ must be.

Each continuity condition, at the knot between two intervals, is a constraint on the B-spline, and the sum of the number of constraints must equal the number of unknown polynomial constants, in order for this B-spline to be well-defined. If the multiplicity of each knot is 1, then the number of constraints is $d$. In general, there are $d+1-m_j$ constraints at each knot, where $m_j$ is knot multiplicity. Over $p$ intervals, that makes $(p-1)(d+1) - \sum_j m_j$ constraints.

There are two other kinds of constraints. We use a normalization of the B-spline to define how tall it will be, which is one constraint. We can also add boundary conditions. Let's first consider the case where the B-spline goes to zero at both sides, with $C^{d-1}$ smoothness and has single multiplicity at each knot. That means there are $d$ boundary conditions at each side. If we balance our unknowns with our constraints, that will tell us how many intervals wide this B-spline extends.

$p(d+1)= (p-1)d + 1 + 2d$

The solution to this is that $p=d+1$. For a smooth polynomial spline, the B-spline for $d=2$ has 3 pieces. For a degree $d=3$ B-spline, there are 4 pieces, and there is a more general form that holds true when there is less smoothness.

We can have up to $2d$ boundary conditions, so let's call that $b\le 2d$. Recall that boundary conditions can be encoded as multiplicity, where a multiplicity of 1 enforces all possible conditions and a multiplicity of $d+1$ enforces no boundary conditions at that side.

$b=2d+2-m_1-m_{p+1}$

We can also have knots with higher multiplicity, which reduces the total number of constraints. Note that the $j$ in this equation runs over internal knots, so it doesn't include endpoints and there will be $(p-1)$ values of $j$.

$p(d+1) = (p-1)(d+1) - \sum_{j=2}^{p} m_j + 1 + b$

If we simplify this equation, it looks like $p$ drops out, but there are $p-1$ values of $m_j$, and usually $m_j=1$ for all knots. The boundary conditions, as well, are $b\le 2d$.

$\sum_{j=2}^{p} m_j = b - d$

As a check, our simple example has multiplicity 1 and $2d=4$ boundary conditions, so $1+1 = 4 - 2$, so math still works. Let's substitute multiplicity at endpoints for $b$.

$\sum_{j=2}^{p} m_j = 2d+2-m_1-m_{p+1} - d$

This simplifies to a formulat that says that, if you represent your boundary conditions and continuity with B-splines, then every consecutive set of points, where the sum of the multiplicity is $d+2$, is the domain of another B-spline.

$\sum_{j=1}^{p+1} m_j = d+2$

We can walk along the array of multiplicity values in order to generate B-spline polynomials. I'll write this as code so that you can see that we walk not by vertex but by multiple of the vertices, as though they were written out as an ordered array that had duplicates.

```julia
function every_bspline_on_axis(x::AbstractArray, multiplicity::AbstractArray)
    # -1 because the last point has no B-spline to the right.
    bspline_count = sum(multiplicity) - 1
    bsplines = []
    vertex_index = 1
    multiplicity_index = 1
    while vertex_index <= bspline_count
        push!(bsplines, generate_bspline(x, multiplicity, vertex_index, multiplicity_index))
        multiplicity_index += 1
        if multiplicity_index > multiplicity[vertex_index]
            vertex_index += 1
            multiplicity_index = 1
        end
    end
end
```

As a mental check, can a single interval have a B-spline? It has no knots. We see that a B-spline, in this case, is well-defined because there is a solution to $0=b-d$, by setting $b$ boundary conditions at the left-hand side $x_l$ and right-hand side $x_r$. For the case of degree 2, there are three ways to do this. 1. Set $f(x_l)=0$, $f'(x_l)=0$. 2. Set $f(x_l)=0$ and $f(x_r)=0$. 3. Set $f(x_r)=0$ and $f'(x_r)=0$. It turns out that *every* interval has $d+1$ splines defined that are nonzero on that interval, no matter what the boundary conditions are.


## How to count B-splines with multiplicity

I'm wondering how to write an algorithm for B-splines that uses multiplicity instead of creating an axis that has repeated values where there are knots. When there is an axis, the answer is simple. If there $k+1$ points on the axis, then there are $k+1-m$ B-splines. The $i$-th B-spline is at axis point $\tau_i$. It's straightforward, but this changes when using multiplicity.

The count of axis points is related to the cumulative multiplicity.

$M_j = \sum_{i=1}^j m_i$

The $i$-th B-spline can be written in terms of multiplicity if we think of it as having $i-1$ cumulative multiplicity to the left of it. If we are looking at interval $j$, then $M_j$ is the total number of axis points to the left. That means the last B-spline to cover this interval is $M_j$ and the first B-spline to cover it will be $M_j-m + 1$, where $m$ is the order.


## Calculating B-splines

We don't calculate the polynomial representation of B-splines. We don't need to, because there are nifty functions that evaluate B-splines without calculating their polynomial constants. The best source for this is Schumaker's *Spline Functions: Basic Theory.*

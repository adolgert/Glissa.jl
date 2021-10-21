# B-splines

## Introduction

### Compare with polynomials

The [`Notation`](@ref) section described polynomial splines as a set of polynomials defined on neighboring intervals. For computational purposes, we don't usually represent polynomial splines as a list of polynomial constants. Instead, we use B-splines, because they are a more compact representation and offer advantages for computations using polynomial splines.

Consider the difference between polynomial splines and B-splines for evaluation, the moment you compute $y$ from $x$. For polynomial splines, there are knots at $x_i$ and polynomial constants $c_{ij}$. In order to compute $y(x)$, you first find the interval $j$ such that $x_j \le x < x_{j+1}$. Then evaluate

$y(x) = \sum_{i=1}^{i=m} c_{ij} x^{i-1}.$

On the other hand for B-splines, each B-spline is some function $B_j(x)$, where $j$ is the index of that B-spline in the interval. Here, you evaluate $y(x)$ with a sum over all B-splines

$y(x) = \sum_j c_j B_j(x).$

Only a few of the $B_j$ contribute to the sum, but it's a simpler evaluation, and there is only one constant for each B-spline. And, as a neat trick, we don't need to store the polynomials that define the $B_j(x)$. They are quick to calculate from the knots.


### Use for interpolation

One use of splines is to interpolate points. Given a set of values $y_k$ at points, $x_k$, what are the constants, $c_j$ for which

$y_k = \sum_j c_j B_j(x_k)?$

This is a set of $k$ equations with $j$ unknowns. If the equations and unknowns match, then this is a matrix equation to solve. For interpolation, finding the $c_j$ from the $(x_k, y_k)$ is considered a prefiltering step for B-splines because it is then the $c_j$ that you save.


## B-spline definition

B-splines are a basis set for polynomial splines on a particular axis, $x_1 < x_2 < x_3 \ldots < x_{k+1}$, with a given set of natural boundary conditions. Because they are a basis set, every allowable polynomial spline with those boundary conditions can be written as a sum of the B-splines for that axis. The set of B-splines sounds useful, so how could we calculate the polynomial representation of any one of the B-splines for a given axis?

If we look for the polynomial pieces that define a single B-spline, we start by picking a left-most point, say $x_i$. Then the B-spline will be zero to the left, positive for some number of intervals, say $p$ intervals. And then it will be zero to the right of that. We can count the unknown values as $p$ intervals of polynomials of degree $d$, so there are $p(d+1)$ unknown values. Oddly, we don't know the number of intervals for this B-spline yet, but continuity conditions will tell us.

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


## Calculating B-splines

We don't calculate the polynomial representation of B-splines. We don't need to, because there are nifty functions that evaluate B-splines without calculating their polynomial constants. Want to derive those here? No, you can suffer through Schumaker's *Spline Functions: Basic Theory.* It's a rite of passage?

In order to test the software, we can calculate the polynomial representation explicitly. We do that by turning the section above into a set of equations and then solving those equations.

First, define the problem. There will be and axis $(x_1,x_2,\ldots x_{k+1})$. Each internal knot of that axis will have a multiplicity, $m_j$. That's the input, and the output will be a set of B-spline polynomials.

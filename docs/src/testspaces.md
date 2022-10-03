# Testing B-splines

We can test well if we focus on theoretical properties that define B-splines. If this set of theoretical properties is mathematically enough to define B-splines, then we don't need any more tests.

## Locality, positivity, and normalization

For both Q-splines and N-splines, they are zero outside of their domain, where the domain for $Q_i^m(x)$ is $\tau_i \le x < \tau_{i+m}$. They must be positive within that domain with normalization for $Q(x)$ such that the first moment is $1/m$.

$\int_{\tau_i}^{\tau_{i+m}}Q_i^m(x)dx = 1/m$

The second moment is also known.

$\int_{\tau_i}^{\tau_{i+m}}xQ_i^m(x)dx = \frac{\sum_{j=i}^{i+m}\tau_j}{m(m+1)}$

and, for the normalized splines, there is a factor that makes the integral smaller for smaller intervals, which can keep numbers in better check. So we can check how the normalized spline relates to the unnormalized.

$N_i^m(x)=(\tau_{i+m}-\tau_i)Q_i^m(x)$

We also know, for normalized splines, that all splines which intersect a point on the axis, will add up to one.

$\sum_{i=j+1-m}^{j}N_i^m(x) = 1$


## Spans a space

The main claim of B-splines is that they are a local basis for a space of polynomial splines. We can test that by constructing random polynomial splines that obey the rules of the axis and knots, then solving for the B-spline constants that fit them. Our test is that the resulting splines fit the values exactly. We need to make polynomial splines first.

### Build a random polynomial spline

These polynomial splines are defined on an axis, $(\tau_1, \tau_2,\ldots,\tau_k)$ where there are a specified number of boundary conditions at each side, all of which, in this case set the value or a derivative to zero. There are also join conditions on the knots, as specified by multiplicity. We can express this as a set of equations. Define a column vector out of the polymonial coefficients for each polynomial piece on an interval, $\vec{c}$. So that's $(c_{11}, c_{21}, c_{31}, c_{12}, c_{22}\ldots)$ for a set of polynomials where interval $j$ is $\sum_i c_i (x-\tau_j)^{i-1}$. Our equations will look like a matrix.

$A\vec{c} = \vec{b}$

The $A$ are continuity equations at boundary conditions and knots. The $\vec{b}$ are all zeros, because these equations assert either that the boundary condition is zero or that two knots are continuous up to some derivative, which relates the knots on the left, again leaving a zero on the right-hand side. OK, we have more unknowns, $\vec{c}$ than we have equations. That makes $A$ an $m\times n$ matrix where $n>m$. How can we generate consistent $c$ values?

The generalized inverse lets us make these. If we take the singular value decomposition, the truncated version with tildes,

$A=\tilde{U}\tilde{\Sigma}\tilde{V}^*$

and take its inverse, we get $A^{\dag}=\tilde{V}\tilde{\Sigma}^{-1}\tilde{U}^*$. Then we can make many $\vec{c}$ by projection.

$\vec{c}=A^{\dag}b + (I-A^{\dag}A)\vec{w}$

The $\vec{w}$ can be any values we want. The first term is a solution to the equations, and the second term is removing from our choice of $\vec{w}$ all contributions that don't solve the equations.

## Find a B-spline approximation to values generated from that polynomial spline

So we can generate a random vector $\vec{w}$ which gives us a polynomial spline defined by $\vec{c}$. At this point, if we can find a linear combination of splines that matches that polynomial spline, then that's a good test. We set up that question as an equation, too.

Here, the unknown $\vec{x}$ is the set of coefficients of all B-splines defined for an axis. For each of these, let's generate a bunch of points from every interval on the axis, If the order is $m$, then generate at least $m$ points for each interval $\tau_i \le x < \tau_{i+1}$. These are our $\vec{y}$ on the right-hand side of an equation.

$M\vec{x}=\vec{y}$

The matrix $M$ is created by making one row for each B-spline and one column for each of the sample points. Then we do a linear regression on this to get a closest set of $\vec{x}$.

## Estimate error between observed and fit

Finally, estimate the error.

$|\vec{y}-M\vec{x}|$

That should be very small for all $\vec{w}$ we can create. We can run this in a loop with randomly-generated polynomial splines.

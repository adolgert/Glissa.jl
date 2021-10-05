# Notation

This package works with polynomials. It will work from traditional polynomial sources, and the notation isn't great. This document will work through what notation to use and how to label variables in the code.

* How to represent points that implicitly extend to either side of the axis. Schumaker talks about an expansion, $y_i$ of the $x_i$. This has been confusing to read.

* What to call each interval - polynomial pieces is simple.

* How to deal with coincident knots. Some algorithms assume consecutive knots can be equal, while others assign a multiplicity to knots. Relying on floating-point equality is a bad idea, so I'd like to use multiplicity.

* There is the order, $m$ and degree $d$. The degree is a number people see when they write out the polynomial, so let's prefer the degree and call it $d$.

The code will need to walk through axis points, taking steps through points that have multiplicity. We can make an inner loop from 1 to m[i] for each point. I think that will work. If we need an index on the virtual points, then that can use an `i+=1` scheme.

# What is a B-spline

## Polynomial splines

Let's start with an axis on the real number line. Set a point $x_1$. A second-degree polynomial is

$c_1 + c_2 (x - x_1) + c_3 (x - x_1)^2.$

This is also called an order 3 polynomial because it has 3 terms. Spline literature likes to specify the order instead of the degree.

If we divide the axis into pieces between some $x_1$ and $x_n$, we get $x_1 < x_2 < x_3 \cdots < x_n$. We can define a separate polynomial piece on each interval. For example, the second piece is a polynomial to the right of $x_2$, $d_1 + d_2 (x - x_2) + d_3 (x - x_2)^2$. Each polynomial piece could be completely separate, not meeting at all. Let's discuss how two polynomials meet.

If the polynomials are continuous across $x_2$ --- the $x$-values where the axis has breaks are called "knots" --- then the two polynomials are related by an equation.

$c_1 + c_2 (x_2 - x_1) + c_3 (x_2 - x_1)^2 = d_1 + d_2 (x_2 - x_2) + d_3 (x_2 - x_2)^2 = d_1$

If the polynomials have the same derivative at $x_2$, then this equation holds.

$c_2 + 2c_3 (x_2 - x_1) = d_2 + 2d_3 (x_2 - x_2)$

If we try to require that the second derivative of a second-degree polynomial matches at $x_2$, then we see that, by this point, the two polynomials have to be the same.

$c_3 = d_3$

That means that a polynomial spline is constructed from pieces which can have no continuity or contiuous derivatives up to a jth derivative, where $j$ is less than the degree. Sometimes we call continuity a 0th derivative continuity.

We're left with polynomial pieces where, at each knot, there is a specified continuity.

## What B-splines are for

Construct a polynomial spline from B-splines. The B-splines are a basis set for polynomial splines. If you had 100 pieces to an interval and wanted to define degree 3 splines, you would need to store 400 polynomial coefficients. But if you use B-splines, you only need to store 103 B-spline coefficients. And, what's more, if you had predetermined that the polynomial spline would have certain continuity conditions between pieces, you'd have to ensure the polynomial coefficients retained that continuity as the values changed, but the B-splines have continuity conditions baked into their selection, so that any set of B-spline coefficients always defines a polynomial spline of the desired continuity.

So we're going to take the degree of the spline as a given, and we're going to take the internal continuity conditions as a given.

## Boundary conditions

At the left-most knot and the right-most knot, we can define boundary conditions. A boundary condition is a specific value for the $j$-th derivative of the polynomial at that boundary. We will discuss so-called natural boundary conditions, meaning that we require the polynomial value to be 0 at the boundary, or its value and first derivative both to be 0, or its value and first two derivatives to be 0, and so on.

## B-spline definition

Given an axis with knots, $x_1 < x_2 < x_3 \cdots < x_n$, every polynomial spline that is uniquely-defined by its continuity and natural boundary conditions is a B-spline. If we make a list of the B-splines, then we can guarantee that any polynomial spine with the same continuity can be described by a sum of B-splines. The B-splines form a basis set for polynomial splines.

A polynomial spline that has no specified internal values is uniquely-defined when the continuity and boundary conditions, together, are enough to specify the constants, $c_1$, $c_2$, of its polynomial pieces. We figure out what these are by playing a game where we a) count the unknowns and b) count the equations that constrain their unknowns. These must be equal in order to define a B-spline.

If the B-spline has degree $d$, then each polynomial piece has $d+1$ unknowns. There will be $p$ polynomial pieces, defined by $p+1$ knots on the axis. The total number of unknowns is therefore $(d+1)p$.

Start with internal continuity. Each continuity level is another constraint. Let's limit ourselves to B-splines that have the same continuity at every internal knot. Where there is maximal continuity, so up-to-the $j$th derivative matches for degree $d=j+1$, we'll call that $m=1$. And when $m=d + 1$, that's no continuity between polynomial pieces. For $p$ pieces, there are $p-1$ places those pieces join internally, which adds $(p-1)(d + 1 - m)$ equations to our constraints.

So far, the score is that there are $(d+1)p$ unknowns and $(p-1)(d + 1 - m)$ constraints, which leaves $d+1+mp-m$ missing constraints, and we get them from two places. One is that there needs to be an overall scaling factor for this polynomial, something to determine its total area. For $Q$ splines, that total area is set at $1/(d+1)$. For normalized splines, that total area is set at $(x_l-x_1)/(d+1)$. We'll see those choices work out nicely. That leaves $d+mp-m$ constraints, which we set as boundary conditions.

There are up to $d$ boundary conditions on each side of the polynomial splines (for degree 2, that's zeroth and first derivatives), so we can define a B-spline as long as

$d + m(p-1) \le 2d.$

As long as that holds, there is some combination of natural boundary conditions that will define a B-spline. If we shift those variables around, we get $m(p-1) \le d$. The most-smooth spline has $m=1$ because it has the most continuity at internal knots, so that $p=d+1$. That means this kind of spline has three polynomial pieces for a second-degree spline and four polynomial pieces for a third-degree spline.

## Use as a basis

For a given spline degree and specified internal continuity, including all b-splines on an axis, including all allowed natural boundary conditions at the edges of the axis, gives a basis for polynomials of that degree and internal continuity.

## Finding B-splines

Using a matrix of equations. Using the magical functions.

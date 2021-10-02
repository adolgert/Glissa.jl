# Notation

This package works with polynomials. It will work from traditional polynomial sources, and the notation isn't great. This document will work through what notation to use and how to label variables in the code.

* How to represent points that implicitly extend to either side of the axis. Schumaker talks about an expansion, $y_i$ of the $x_i$. This has been confusing to read.

* What to call each interval - polynomial pieces is simple.

* How to deal with coincident knots. Some algorithms assume consecutive knots can be equal, while others assign a multiplicity to knots. Relying on floating-point equality is a bad idea, so I'd like to use multiplicity.

* There is the order, $m$ and degree $d$. The degree is a number people see when they write out the polynomial, so let's prefer the degree and call it $d$.

The code will need to walk through axis points, taking steps through points that have multiplicity. We can make an inner loop from 1 to m[i] for each point. I think that will work. If we need an index on the virtual points, then that can use an `i+=1` scheme.


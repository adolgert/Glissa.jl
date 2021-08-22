# Glissa

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://adolgert.github.io/Glissa.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://adolgert.github.io/Glissa.jl/dev)
[![Build Status](https://github.com/adolgert/Glissa.jl/workflows/CI/badge.svg)](https://github.com/adolgert/Glissa.jl/actions)
[![Coverage](https://codecov.io/gh/adolgert/Glissa.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/adolgert/Glissa.jl)

This has a Hermite cubic spline interpolator. It does monotonic fitting using
Hyman splines.

```{Julia}
x = [7.99, 8.09, 8.19, 8.7, 9.2, 10, 12, 15, 20]
f = [0, 2.76429e-5, 4.37498e-2, 0.169183, 0.469428, 0.943740, 0.998636, 0.999919, 0.999994]
endpoint_derivatives = [0.0, 0.0]
cs = cubic_spline(x, f, endpoint_derivatives)
y = cs.(8:0.01:20)
```

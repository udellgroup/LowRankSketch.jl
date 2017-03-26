# LowRankSketch

[![Build Status](https://travis-ci.org/madeleineudell/LowRankSketch.jl.svg?branch=master)](https://travis-ci.org/madeleineudell/LowRankSketch.jl)

[![Coverage Status](https://coveralls.io/repos/madeleineudell/LowRankSketch.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/madeleineudell/LowRankSketch.jl?branch=master)

[![codecov.io](http://codecov.io/github/madeleineudell/LowRankSketch.jl/coverage.svg?branch=master)](http://codecov.io/github/madeleineudell/LowRankSketch.jl?branch=master)

This package implements a forthcoming paper by Tropp, Yurtsever, Udell and Cevher on sketching methods for low rank matrix approximation.

# Syntax

This code uses standard Julia v0.6 syntax

* a `type` is a container for data
* `<:` and `::` are type assertions; read them as "is a"
* a `type` can inherit from an `abstract type`. For example, the type statement
```
abstract type Bar end
type Foo<:bar end
```
means that the `type` `Foo` inherits from the `type` `Bar`
* If methods are not defined for the `type`, then they fall back to the default methods defined for the `abstract type`. For example,
```
function say_hello(thing::Bar)
  print("hello, bar")
end
f = Foo()
say_hello(f)
```
will print `hello, bar`: `f::Bar` is true since `Foo` inherits from `Bar`.
* curly brackets `{}` after a type or function name tell you that the types in curly brackets parametrize the type or function. So, eg, `Matrix{Float64}` is a matrix of floating point numbers with 64 bit precision

# Differences with pseudocode

* To avoid defining ctranspose for a DimRedux map, I've moved the ctranspose inside right-multiplication by a DimRedux map.
* I think there's something wrong in the definition of multiplication for SSFRT; the dimensions don't even work.

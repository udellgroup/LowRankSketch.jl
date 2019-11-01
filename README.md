# LowRankSketch

[![Build Status](https://travis-ci.org/madeleineudell/LowRankSketch.jl.svg?branch=master)](https://travis-ci.org/madeleineudell/LowRankSketch.jl)

[![Coverage Status](https://coveralls.io/repos/madeleineudell/LowRankSketch.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/madeleineudell/LowRankSketch.jl?branch=master)

[![codecov.io](http://codecov.io/github/madeleineudell/LowRankSketch.jl/coverage.svg?branch=master)](http://codecov.io/github/madeleineudell/LowRankSketch.jl?branch=master)

This package implements a forthcoming paper by Tropp, Yurtsever, Udell and Cevher on sketching methods for low rank matrix approximation.

# Syntax

This code uses standard Julia v0.6 syntax

* a `type` is a container for data
* `<:` and `::` are type assertions; read them as "is a"
* a type can inherit from an abstract type. For example, the type statement
```julia
abstract type Bar end
type Foo<:Bar end
```
means that the type `Foo` inherits from the type `Bar`.
* functions are defined to take arguments of specific types. So, if we define
```julia
function say_hi(thing::Foo)
  print("hi, foo")
end
```
and call `say_hi(Foo())`, it will correctly print `"hi, foo"`.
But `say_hi(5)` will be an error, since `5` is not a `Foo`.
* If methods are not defined for the type, then they fall back to the default methods defined for the abstract type. For example,
```julia
function say_hello(thing::Bar)
  print("hello, bar")
end
f = Foo()
say_hello(f)
```
will print `hello, bar`. (`f::Bar` is true since `Foo` inherits from `Bar`.)
* curly brackets `{}` after a type or function name tell you that the types in curly brackets parametrize the type or function. So, eg, `Matrix{Float64}` is a matrix of floating point numbers with 64 bit precision. These are often used in this code to ensure that the field F (say, real or complex numbers) over which the (many) arguments of the function are defined is the same.

# Differences with pseudocode

* The SparseDimRedux is not yet implemented
* The PSDSketch is not yet implemented

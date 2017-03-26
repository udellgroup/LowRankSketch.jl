## Tests for dimension reductions
using LowRankSketch
using FactCheck

TOL = 1e-1

m,n,k,s = 100,200,10,20
m,n,k,s = 10,20,3,3

facts("DimRedux") do

  # real valued Gaussian dimension reduction
  dr = GaussDimRedux(k,n)
  M = rand(n,m)
  dr*M
  M'*dr

  # complex valued Gaussian dimension reduction
  dr = GaussDimRedux(k,n,Complex)
  M = rand(n,m) + im*rand(n,m)
  dr*M
  M'*dr

  # real valued sparse dimension reduction
  dr = SparseDimRedux(k,n)
  M = rand(n,m)
  dr*M
  M'*dr

  # complex valued sparse dimension reduction
  dr = SparseDimRedux(k,n,.2,Complex)
  M = rand(n,m) + im*rand(n,m)
  dr*M
  M'*dr

  # real valued SSRFT dimension reduction
  dr = SSRFTDimRedux(k,n)
  M = rand(n,m)
  dr*M
  M'*dr

  # complex valued SSFRT dimension reduction
  dr = SSRFTDimRedux(k,n,Complex)
  M = rand(n,m) + im*rand(n,m)
  dr*M
  M'*dr

end

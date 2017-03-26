## Tests for correct sketching
using LowRankSketch
using FactCheck

TOL = 1e-1

# we'll sketch rank-1 matrices with a two dimensional sketch
m,n,k,s = 100,200,2,5
# m,n,k,s = 10,20,3,3
srand(1)

facts("Sketch") do

  context("sketch real matrix with GaussDimRedux") do
    sketch = Sketch(m,n,k,s,GaussDimRedux,Float64)
    M = rand(m,1)*rand(1,n)
    linear_update(sketch, M)
    Q,W,P = low_rank_approx(sketch)
    @fact vecnorm(M - Q*W*P')/vecnorm(M) --> roughly(0, TOL)
    U,Σ,V = fixed_rank_approx(sketch,1)
    @fact vecnorm(M - U*Σ*V')/vecnorm(M) --> roughly(0, TOL)
  end

  context("sketch real matrix with SSRFTDimRedux") do
    sketch = Sketch(m,n,k,s,SSRFTDimRedux,Float64)
    M = rand(m,1)*rand(1,n)
    linear_update(sketch, M)
    Q,W,P = low_rank_approx(sketch)
    @fact vecnorm(M - Q*W*P')/vecnorm(M) --> roughly(0, TOL)
    U,Σ,V = fixed_rank_approx(sketch,1)
    @fact vecnorm(M - U*Σ*V')/vecnorm(M) --> roughly(0, TOL)
  end

  context("sketch complex matrix with GaussDimRedux") do
    sketch = Sketch(m,n,k,s,GaussDimRedux,Complex{Float64})
    M = rand(m,1)*rand(1,n) + im*rand(m,1)*rand(1,n)
    linear_update(sketch, M)
    Q,W,P = low_rank_approx(sketch)
    @fact vecnorm(M - Q*W*P')/vecnorm(M) --> roughly(0, TOL)
    U,Σ,V = fixed_rank_approx(sketch,1)
    @fact vecnorm(M - U*Σ*V')/vecnorm(M) --> roughly(0, TOL)
  end

  context("sketch complex matrix with SSRFTDimRedux") do
    sketch = Sketch(m,n,k,s,SSRFTDimRedux,Complex{Float64})
    M = rand(m,1)*rand(1,n) + im*rand(m,1)*rand(1,n)
    linear_update(sketch, M)
    Q,W,P = low_rank_approx(sketch)
    @fact vecnorm(M - Q*W*P')/vecnorm(M) --> roughly(0, TOL)
    U,Σ,V = fixed_rank_approx(sketch,1)
    @fact vecnorm(M - U*Σ*V')/vecnorm(M) --> roughly(0, TOL)
  end

  context("sketch integer matrix with SSRFTDimRedux") do
    # (make sure to use sketch compatible with type of matrix to be sketched, eg integer matrix -> use SSRFT)
    # nb linear update doesn't work yet for this (InexactError) b/c output of SSRFT is not integral
    sketch = Sketch(m,n,k,s,SSRFTDimRedux,Complex{Int})
    sketch = Sketch(m,n,k,s,SSRFTDimRedux,Int)
  end

end

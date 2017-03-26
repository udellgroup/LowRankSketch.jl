export DimRedux, GaussDimRedux, SparseDimRedux, SSRFTDimRedux, *

### Dimension reduction operators

"""DimRedux is an abstract dimension reduction map
over the field F, which is either the real or complex numbers"""
abstract type DimRedux{F<:Number} end

"""By default, a DimRedux map reduces a matrix of data via left or right multiplication by its internal variable Ξ.
We extend the multiplication method * to allow multiplication by dimension reduction maps"""
*{F<:Number}(dr::DimRedux{F}, M::Matrix{F}) = dr.Ξ * M
*{F<:Number}(M::Matrix{F}, dr::DimRedux{F}) = M * dr.Ξ' #' nb: different from paper

"""GaussDimRedux is a kind of dimension reduction map
It reduces data over the field F using the Gaussian random matrix Ξ"""
type GaussDimRedux{F<:Number}<:DimRedux{F}
    Ξ::Matrix{F}
end

"""Constructs a GaussDimRedux map which reduces data in F^n to F^k
By default, the field F is a real (floating point) number"""
function GaussDimRedux{F<:Number}(k::Int, n::Int, ::Type{F}=Float64)
    if F<:Real
        Ξ = randn(k,n)
    elseif F<:Complex
        Ξ = randn(k,n) + im*randn(k,n)
    else error("GaussDimRedux not defined over the field $F")
    end
    return GaussDimRedux(Ξ)
end

"""SparseDimRedux is a kind of dimension reduction map
It reduces data over the field F using the sparse random matrix Ξ"""
type SparseDimRedux{F<:Number}<:DimRedux{F}
    Ξ::SparseMatrixCSC{F}
end

"""Constructs a SparseDimRedux map which reduces data in F^n to F^k
By default, the field F is the real numbers
p is the sparsity level, with default level .1 (???)"""
function SparseDimRedux{F<:Number}(k::Int, n::Int, p::Float64=.1, ::Type{F}=Real)
    if F<:Real
        Ξ = sprandn(k,n,p) # nb --- not defined in paper
    elseif F<:Complex
        Ξ = sprandn(k,n,p) + im*sprandn(k,n,p) # nb --- not defined in paper
    else error("SparseDimRedux not defined over the field $F")
    end
    return SparseDimRedux(Ξ)
end

"""SSRFTDimRedux is a kind of dimension reduction map
It reduces data over the field F using a scrambled SRFT
Its internal variables are
the permutations Π1 and Π2, the random signs ϵ1 and ϵ2, and the coordinate restriction R"""
type SSRFTDimRedux{F<:Number}<:DimRedux{F}
    Π1::Array{Int,1} # first permutation
    ϵ1::Array{F,1} # first random signs
    Π2::Array{Int,1} # second permutation
    ϵ2::Array{F,1} # second random signs
    R::Array{Int,1} # restriction onto k coordinates
end

"""Constructs a SSRFTDimRedux map which reduces data in F^n to F^k
By default, the field F is the real numbers"""
function SSRFTDimRedux{F<:Number}(k::Int, n::Int, ::Type{F}=Float64)
    Π1 = randperm(n) # a random permutation of {1,...,n}
    Π2 = randperm(n) # a random permutation of {1,...,n}
    R = randperm(n)[1:k] # a random selection of k elements from {1,...,n}
    if F<:Real
        ϵ1 = 1.0*one(F) * rand([-1, 1], n) # a random ±1 vector of length n in the field F
        ϵ2 = 1.0*one(F) * rand([-1, 1], n) # a random ±1 vector of length n in the field F
    elseif F<:Complex
        ϵ1 = 1.0*one(F) * rand([-1, 1], n) + im * 1.0 * one(F) * rand([-1, 1], n) # a random ±1±i vector of length n in the field F
        ϵ2 = 1.0*one(F) * rand([-1, 1], n) + im * 1.0 * one(F) * rand([-1, 1], n) # a random ±1±i vector of length n in the field F
    else error("SSRFTDimRedux not defined over the field $F")
    end
    return SSRFTDimRedux(Π1,ϵ1,Π2,ϵ2,R)
end

"""The SSRFTDimRedux map uses an efficient specialized method to reduce a matrix of data,
multiplying it on the left by Ξ = R FFT ϵ2 Π2 FFT ϵ1 Π1"""
function *{F<:Number}(dr::SSRFTDimRedux{F}, M::Matrix{F})
    B = M[dr.Π1,:] # scramble rows
    B .= dr.ϵ1 .* B # scramble signs - equivalent to diagm(dr.ϵ1) * B
    if F<:Real
        B = dct(B)
    elseif F<:Complex
        B = fft(B)
    end
    B = B[dr.Π2,:] # scramble rows
    B .= dr.ϵ2 .* B # scramble signs - equivalent to diagm(dr.ϵ2) * B
    if F<:Real
        B = dct(B)
    elseif F<:Complex
        B = fft(B)
    end
    B = B[dr.R,:] # select rows
    return B
end
"""The SSRFTDimRedux map uses an efficient specialized method to reduce a matrix of data,
multiplying it on the right by Ξ' = (R FFT ϵ2 Π2 FFT ϵ1 Π1)'"""
function *{F<:Number}(M::Matrix{F}, dr::SSRFTDimRedux{F})
  B = M[:,dr.Π1] # scramble columns
  B .= B .* dr.ϵ1' # scramble signs - equivalent to B * diagm(dr.ϵ1)
  if F<:Real
      B = dct(B)
  elseif F<:Complex
      B = ifft(B)
  end
  B = M[:,dr.Π2] # scramble columns
  B .= B .* dr.ϵ2' # scramble signs - equivalent to B * diagm(dr.ϵ1)
  if F<:Real
      B = dct(B)
  elseif F<:Complex
      B = ifft(B)
  end
  B = B[:,dr.R] # select columns
  return B
end

# """Constructs a SSRFTDimRedux map which reduces data in F^n to F^k
# By default, the field F is the real numbers"""
# function SSRFTDimRedux{F<:Number}(k::Int, n::Int, ::Type{F}=Real)
#     perm = randperm(n) # a random permutation of {1,...,n}
#     coords = randperm(n)[1:k] #  k elements chosen without replacement from {1,...,n}
#     if F<:Real
#         ϵ = one(F) * rand([-1, 1], n) # a random ±1 vector of length n in the field F
#     elseif F<:Complex
#         ϵ = one(F) * rand([-1, 1], n) + im * one(F) * rand([-1, 1], n) # a random ±1±i vector of length n in the field F
#     else error("SSRFTDimRedux not defined over the field $F")
#     end
#     return SSRFTDimRedux(perm,ϵ,coords)
# end
#
# """The SSRFTDimRedux map uses an efficient specialized method to reduce a matrix of data,
# rather than the default method of left and right multiplication by Ξ"""
# function *{F<:Number}(dr::SSRFTDimRedux{F}, M::Matrix{F})
#     B = M[dr.perm,:] # scramble rows
#     B = diagm(dr.ϵ) * B # scramble signs
#     if F<:Real
#         B = dct(B)
#     elseif F<:Complex
#         B = fft(B)
#     end
#     B = B[dr.coords,:] # select rows --- nb: different from paper
#     return B
# end
# function *{F<:Number}(M::Matrix{F}, dr::SSRFTDimRedux{F})
#     B = M[:,dr.perm] # scramble columns
#     B = B * diagm(dr.ϵ) # scramble signs --- nb: different from paper
#     if F==Real
#         B = dct(B)
#     elseif F==Complex
#         B = fft(B)
#     end
#     B = B[:,dr.coords] # select columns --- nb: different from paper
#     return B
# end

export Sketch, linear_update, low_rank_approx, fixed_rank_approx

### Sketching operators and methods

"""Sketch is a sketch for compressing a matrix over the field F
It is parametrized by the field F, which may be either the real or complex numbers,
and by the type of dimension reduction map DR, which may be either GaussDimRedux, SparseDimRedux, or SSRFTDimRedux
The matrices X, Y and Z and the dimension reduction map are defined over the same field F"""
type Sketch{F<:Number,DR<:DimRedux{F}}
    Υ::DR        # from F^m to F^k
    Ω::DR        # from F^n to F^k
    Φ::DR        # from F^m to F^s
    Ψ::DR        # from F^n to F^s
    X::Matrix{F}  # in F^{k×n}
    Y::Matrix{F}  # in F^{m×k}
    Z::Matrix{F}  # in F^{s×s}
end

"""Constructs a sketch for a matrix in F^{m×n}
using the dimension reduction map DR.
By default, the field F is the real (floating point) numbers"""
function Sketch{F<:Number,DR<:DimRedux}(m::Int, n::Int, k::Int, s::Int, ::Type{DR}, ::Type{F}=Float64) # DR<:DimRedux{F}??
    Υ = DR(k,m,F)
    Ω = DR(k,n,F)
    Φ = DR(s,m,F)
    Ψ = DR(s,n,F)
    X = zeros(F,k,n)
    Y = zeros(F,m,k)
    Z = zeros(F,s,s)
    return Sketch(Υ,Ω,Φ,Ψ,X,Y,Z)
end

function linear_update{F<:Number,DR<:DimRedux{F}}(sk::Sketch{F,DR}, H::Matrix{F}, θ::F=one(F), τ::F=one(F))
    # the dots are julia syntax to improve memory efficiency
    # nb removed transposes here; put them in the definition of left/right multiplication
    sk.X .= θ*sk.X .+ τ*(sk.Υ*H)
    sk.Y .= θ*sk.Y .+ τ*(H*sk.Ω)
    sk.Z .= θ*sk.Z .+ τ*(sk.Φ*H*sk.Ψ)
    return sk
end

function low_rank_approx(sk::Sketch)
    P,_ = qr(sk.X')
    Q,_ = qr(sk.Y)
    U1,T1 = qr(sk.Φ*Q)
    U2,T2 = qr(sk.Ψ*P)
    # these work for real gaussian case
    # W = pinv(T1)*(U1'*sk.Z*U2)*pinv(conj(T2'))
    # W = T1\(U1'*sk.Z*U2)/(T2')
    W = T1\(U1'*sk.Z*U2)/(T2')
    return (Q,W,P)
end

function fixed_rank_approx(sk::Sketch, r::Int)
    Q,W,P = low_rank_approx(sk)
    svdW,_ = svds(W,nsv=r)
    U,Σ,V = svdW.U, svdW.S, svdW.Vt'
    return (Q*U,diagm(Σ),P*V)
end

import Base.show
export show

function show{F}(io::IO, dr::DimRedux{F})
    k,n = size(dr.Ξ)
    print(io, "$(typeof(dr)): dimension reduction map over the field $F from $n to $k dimensions")
end

function show{F,DR}(io::IO, sk::Sketch{F,DR})
    k,n = size(sk.X)
    m,k = size(sk.Y)
    s,s = size(sk.Z)
    print(io, "$(typeof(sk)): sketch for an $m×$n matrix over the field $F using a $DR dimension reduction map")
end

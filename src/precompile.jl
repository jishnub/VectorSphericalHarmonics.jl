precompile(genspharm, (Int, Int, Float64, Float64))
for LT in [ZeroTo{false}, UnitRange{Int}], MT in [ZeroTo{true}, FullRange{true}, UnitRange{Int}]
    precompile(VSHCache, (Type{Float64}, ML{LT, MT}))
    precompile(VSHCache, (Type{Float64}, LM{LT, MT}))
    for YT in [PB, Hansen, Irreducible], B in [Cartesian, Polar, SphericalCovariant, HelicityCovariant]
        precompile(VSHCache, (Type{Float64}, YT, B, ML{LT, MT}))
        precompile(VSHCache, (Type{Float64}, YT, B, LM{LT, MT}))
    end
end

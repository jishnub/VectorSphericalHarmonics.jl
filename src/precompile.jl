precompile(genspharm, (Int, Int, Float64, Float64))
for LT in [ZeroTo{false}, UnitRange{Int}], MT in [ZeroTo{true}, FullRange{true}, UnitRange{Int}]
    precompile(VSHCache, (Type{Float64}, ML{LT, MT}))
    precompile(VSHCache, (Type{Float64}, LM{LT, MT}))
    lr = LT === UnitRange{Int} ? LT(0:0) : LT(0)
    mr = MT === UnitRange{Int} ? MT(0:0) : MT(0)
    C = VectorSphericalHarmonics.VSHCache(Float64, LM(lr, mr))
    precompile(genspharm!, (typeof(C), Float64, Float64))
    C = VectorSphericalHarmonics.VSHCache(Float64, ML(lr, mr))
    precompile(genspharm!, (typeof(C), Float64, Float64))
    for YT in [PB, Hansen, Irreducible], B in [Cartesian, Polar, SphericalCovariant, HelicityCovariant]
        precompile(VSHCache, (Type{Float64}, YT, B, ML{LT, MT}))
        precompile(VSHCache, (Type{Float64}, YT, B, LM{LT, MT}))
        C = VectorSphericalHarmonics.VSHCache(Float64, YT(), B(), LM(lr, mr))
        precompile(vshbasis!, (typeof(C), YT, B, Float64, Float64))
        C = VectorSphericalHarmonics.VSHCache(Float64, YT(), B(), ML(lr, mr))
        precompile(vshbasis!, (typeof(C), YT, B, Float64, Float64))
    end
end

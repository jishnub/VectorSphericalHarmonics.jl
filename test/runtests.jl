using VectorSphericalHarmonics
using Test
using Aqua
using LinearAlgebra
using LegendrePolynomials
using OffsetArrays
using SphericalHarmonics
using SphericalHarmonics: NorthPole, SouthPole, getY, eltypeY
using SphericalHarmonicModes
using WignerD
using HCubature
using VectorSphericalHarmonics: basisconversionmatrix, cache, cache!
using Rotations
using StaticArrays

@testset "project quality" begin
    if VERSION >= v"1.6"
        Aqua.test_all(VectorSphericalHarmonics, ambiguities = (recursive = false,))
    else
        Aqua.test_all(VectorSphericalHarmonics, ambiguities = false)
    end
end

isapproxdefault(a, b) = isapprox(a, b, atol = 1e-14, rtol = sqrt(eps(Float64)))
isapproxdefault(x::Tuple{Any, Any}) = isapproxdefault(x...)

function genspharm2(l, m, n, θ, ϕ)
    @assert abs(m) <= l "m must satisfy -l <= m <= l, received l = $l, m = $m"
    @assert abs(n) <= l "n must satisfy -l <= n <= l, received l = $l, n = $n"
    √((2l+1)/4π) * WignerD.wignerdjmn(l, m, n, θ) * cis(m * ϕ)
end

const S = cache(10);
const Ylm = getY(S);

@testset "Basis" begin
    for θ in LinRange(0, pi, 10), ϕ in LinRange(0, 2pi, 10)
        for B2 in [HelicityCovariant(), Polar(), SphericalCovariant(), Cartesian()]
            for B1 in [HelicityCovariant(), SphericalCovariant()]
                M = basisconversionmatrix(B1, B2, θ, ϕ)
                Minv = basisconversionmatrix(B2, B1, θ, ϕ)
                # transformations must be unitary
                @test M' * M ≈ I
                @test Minv ≈ M'
                @test abs(det(M)) ≈ 1
            end
        end
    end
    @testset "multiple transformations" begin
        for θ in LinRange(0, pi, 10), ϕ in LinRange(0, 2pi, 10)
            for B in [HelicityCovariant(), Polar(), SphericalCovariant(), Cartesian()]
                M1 = basisconversionmatrix(HelicityCovariant(), B, θ, ϕ)
                M2 = basisconversionmatrix(B, SphericalCovariant(), θ, ϕ)
                M3 = basisconversionmatrix(HelicityCovariant(), SphericalCovariant(), θ, ϕ)
                @test begin
                    res = M2 * M1 ≈ M3
                    if !res
                        println(HelicityCovariant(),"->",B,"->",SphericalCovariant(),", θ = $θ, ϕ = $ϕ")
                    end
                    res
                end

                M1 = basisconversionmatrix(SphericalCovariant(), B, θ, ϕ)
                M2 = basisconversionmatrix(B, HelicityCovariant(), θ, ϕ)
                M3 = basisconversionmatrix(SphericalCovariant(), HelicityCovariant(), θ, ϕ)
                @test begin
                    res = M2 * M1 ≈ M3
                    if !res
                        println(SphericalCovariant(),"->",B,"->",HelicityCovariant(),", θ = $θ, ϕ = $ϕ")
                    end
                    res
                end
            end
        end
    end
    @testset "basis at poles" begin
        @testset "Cartesian and polar" begin
            M = basisconversionmatrix(Cartesian(), Polar(), 0, 0)
            @test M * [1,0,0] ≈ [0,1,0] # x̂ == θ̂ at the north pole
            @test M * [0,1,0] ≈ [0,0,1] # ŷ == ϕ̂ at the north pole
            @test M * [0,0,1] ≈ [1,0,0] # ẑ == r̂ at the north pole

            M = basisconversionmatrix(Polar(), Cartesian(), pi, 0)
            @test M * [1,0,0] ≈ [0,0,-1] # r̂ == -ẑ at the south pole
            @test M * [0,1,0] ≈ [-1,0,0] # θ̂ == -x̂ at the south pole
            @test M * [0,0,1] ≈ [0,1,0]  # ϕ̂ == ŷ at the south pole
        end
    end
end

@testset "cache" begin
    θ, ϕ = pi/3, pi/4
    S2 = @inferred cache(Float64, 1)
    cache!(S2, θ, ϕ)
    @test getY(S2) == getY(cache(Float64, θ, ϕ, 1))

    S2 = @inferred cache(1)
    cache!(S2, θ, ϕ)
    @test getY(S2) == getY(cache(Float64, θ, ϕ, 1))

    @test vshbasis(PB(), Polar(), 1, 1, θ, ϕ, S2) == vshbasis(PB(), Polar(), 1, 1, θ, ϕ)
    @test genspharm(1, 1, θ, ϕ, S2) == genspharm(1, 1, θ, ϕ)
end

@testset "VSH orientations" begin
    θ1, ϕ1 = pi/3, pi/3;
    θNP, ϕNP = 0, 0;
    @testset "PB" begin
        @testset "HelicityCovariant" begin
            cache!(S, θ1, ϕ1, 0)
            Y = @inferred vshbasis(PB(), HelicityCovariant(), 0, 0, θ1, ϕ1, S)
            # for l=0 there is only one non-zero vector Ylm0, which is radial
            for n in axes(Y,2), basisind in axes(Y,1)
                if !(n == 0 && basisind == 0 #= e₀′ (radial) component =#)
                    @test isapproxdefault(0, Y[basisind, n])
                end
            end
        end
        @testset "Polar" begin
            cache!(S, θ1, ϕ1, 0)
            Y = @inferred vshbasis(PB(), Polar(), 0, 0, θ1, ϕ1, S)
            # for l=0 there is only one non-zero vector Ylm0, which is radial
            for n in axes(Y,2), basisind in axes(Y,1)
                if !(n == 0 && basisind == 1 #= radial component =#)
                    @test isapproxdefault(0, Y[basisind, n])
                end
            end
        end
        @testset "SphericalCovariant" begin
            # north pole
            cache!(S, θNP, ϕNP, 0)
            Y = @inferred vshbasis(PB(), SphericalCovariant(), 0, 0, θNP, ϕNP, S)
            # for l=0 there is only one non-zero vector Ylm0, which is radial
            # at north pole r̂ == ẑ. so there is only one non-zero component of Ylm0
            for n in axes(Y,2), basisind in axes(Y,1)
                if !(n == 0 && basisind == 0 #= e₀ (z) component =#)
                    @test isapproxdefault(0, Y[basisind, n])
                end
            end
        end
        @testset "Cartesian" begin
            # north pole
            cache!(S, θNP, ϕNP, 0)
            Y = @inferred vshbasis(PB(), Cartesian(), 0, 0, θNP, ϕNP, S)
            # for l=0 there is only one non-zero vector Ylm0, which is radial
            # at north pole r̂ == ẑ. so there is only one non-zero component of Ylm0
            for n in axes(Y,2), basisind in axes(Y,1)
                if !(n == 0 && basisind == 3 #= z component =#)
                    @test isapproxdefault(0, Y[basisind, n])
                end
            end
        end
    end
    @testset "Hansen" begin
        @testset "HelicityCovariant" begin
            cache!(S, θ1, ϕ1, 0)
            Y = @inferred vshbasis(Hansen(), HelicityCovariant(), 0, 0, θ1, ϕ1, S)
            # for l=0 there is only one non-zero vector Ylm⁽⁻¹⁾, which is radial
            for n in axes(Y,2), basisind in axes(Y,1)
                if !(n == -1 && basisind == 0 #= e₀′ (radial) component =#)
                    @test isapproxdefault(0, Y[basisind, n])
                end
            end
        end
        @testset "Polar" begin
            cache!(S, θ1, ϕ1, 0)
            Y = @inferred vshbasis(Hansen(), Polar(), 0, 0, θ1, ϕ1, S)
            # for l=0 there is only one non-zero vector Ylm⁽⁻¹⁾, which is radial
            for n in axes(Y,2), basisind in axes(Y,1)
                if !(n == -1 && basisind == 1 #= radial component =#)
                    @test isapproxdefault(0, Y[basisind, n])
                end
            end
        end
        @testset "SphericalCovariant" begin
            # north pole
            cache!(S, θNP, ϕNP, 0)
            Y = @inferred vshbasis(Hansen(), SphericalCovariant(), 0, 0, θNP, ϕNP, S)
            # for l=0 there is only one non-zero vector Ylm⁽⁻¹⁾, which is radial
            # at north pole r̂ == ẑ. so there is only one non-zero component of Ylm0
            for n in axes(Y,2), basisind in axes(Y,1)
                if !(n == -1 && basisind == 0 #= e₀ (z) component =#)
                    @test isapproxdefault(0, Y[basisind, n])
                end
            end
        end
        @testset "Cartesian" begin
            # north pole
            cache!(S, θNP, ϕNP, 0)
            Y = @inferred vshbasis(Hansen(), Cartesian(), 0, 0, θNP, ϕNP, S)
            # for l=0 there is only one non-zero vector Ylm⁽⁻¹⁾⁾, which is radial
            # at north pole r̂ == ẑ. so there is only one non-zero component of Ylm0
            for n in axes(Y,2), basisind in axes(Y,1)
                if !(n == -1 && basisind == 3 #= z component =#)
                    @test isapproxdefault(0, Y[basisind, n])
                end
            end
        end
    end
    @testset "Irreducible" begin
        cache!(S, θ1, ϕ1, 0)
        Y = @inferred vshbasis(Irreducible(), Polar(), 0, 0, θ1, ϕ1, S)
        @test all(iszero, Y[:, -1:0])
        @test all(x -> isapproxdefault(x, 0), Y[2:3, 1])
        @test Y[1, 1] ≈ -1/√(4pi)

        for m = -1:1
            Y = @inferred vshbasis(Irreducible(), SphericalCovariant(), 1, m, θ1, ϕ1, S)
            @test Y[m, 0] ≈ 1/√(4pi)
            @test all(iszero, Y[-1:m-1, 0])
            @test all(iszero, Y[m+1:1, 0])
        end

        cache!(S, θNP, ϕNP, 0)
        Y = @inferred vshbasis(Irreducible(), Cartesian(), 0, 0, θNP, ϕNP, S)
        @test all(iszero, Y[:, -1:0])
        @test all(iszero, Y[1:2, 1])
        @test Y[3, 1] ≈ -1/√(4pi)
    end
end

@testset "norm in different bases" begin
    lmax = 4
    for θ in LinRange(0, pi, 10), ϕ in LinRange(0, 2pi, 10)
        cache!(S, θ, ϕ)
        for l = 0:lmax, m = -l:l
            Y = vshbasis(Irreducible(), SphericalCovariant(), l, m, θ, ϕ, S)
            Yp = parent(Y)
            for B in [HelicityCovariant(), Polar(), Cartesian()]
                Y2 = vshbasis(Irreducible(), B, l, m, θ, ϕ, S)
                Y2p = parent(Y2)
                for l in 1:3
                    @test norm(Y2p[:, l]) ≈ norm(Yp[:, l])
                end
            end
        end
    end
end

@testset "VSH components" begin
    lmax = 5
    @testset "HelicityCovariant" begin
        @testset "PB" begin
            for θ in LinRange(0, pi, 10), ϕ in LinRange(0, 2pi, 10)
                cache!(S, θ, ϕ)
                for l in 0:lmax, m in -l:l
                    Y = vshbasis(PB(), HelicityCovariant(), l, m, θ, ϕ, S)
                    @test isapproxdefault(Y[0,0], Ylm[(l,m)])
                    for n in -1:1
                        if n in -min(1, l):min(1, l)
                            @test isapproxdefault(Y[n,n], genspharm2(l, m, n, θ, ϕ))
                        else
                            @test isapproxdefault(Y[n,n], 0)
                        end
                        for basisind in -1:1
                            basisind == n && continue
                            @test isapproxdefault(Y[basisind,n], 0)
                        end
                    end
                end
            end
        end
        @testset "Hansen" begin
            # l = 0 has only one component
            for θ in LinRange(0, pi, 10), ϕ in LinRange(0, 2pi, 10)
                cache!(S, θ, ϕ, 0)
                Y = vshbasis(Hansen(), HelicityCovariant(), 0, 0, θ, ϕ, S)
                @test isapproxdefault(Y[ 0, -1],  √(1/4pi) * WignerD.wignerDjmn(0, 0, 0, 0, θ, ϕ))
                for n in -1:1, basisind in -1:1
                    basisind == 0 && n == -1 && continue
                    @test isapproxdefault(Y[basisind, n], 0)
                end
            end
            for θ in LinRange(0, pi, 10), ϕ in LinRange(0, 2pi, 10)
                cache!(S, θ, ϕ)
                for l in 1:lmax, m in -l:l
                    norm = √((2l+1)/8pi)
                    sqrt2norm = √((2l+1)/4pi)
                    Y = vshbasis(Hansen(), HelicityCovariant(), l, m, θ, ϕ, S)
                    @test isapproxdefault(Y[-1, +1],  norm * WignerD.wignerDjmn(l,  1, -m, 0, θ, ϕ))
                    @test isapproxdefault(Y[ 0, +1],  0)
                    @test isapproxdefault(Y[+1, +1],  norm * WignerD.wignerDjmn(l, -1, -m, 0, θ, ϕ))
                    @test isapproxdefault(Y[-1,  0],  norm * WignerD.wignerDjmn(l,  1, -m, 0, θ, ϕ))
                    @test isapproxdefault(Y[ 0,  0],  0)
                    @test isapproxdefault(Y[+1,  0], -norm * WignerD.wignerDjmn(l, -1, -m, 0, θ, ϕ))
                    @test isapproxdefault(Y[-1, -1],  0)
                    @test isapproxdefault(Y[ 0, -1],  sqrt2norm * WignerD.wignerDjmn(l, 0, -m, 0, θ, ϕ))
                    @test isapproxdefault(Y[+1, -1],  0)
                end
            end
        end
        @testset "Irreducible" begin
            # j = 0 has only one component
            for θ in LinRange(0, pi, 10), ϕ in LinRange(0, 2pi, 10)
                cache!(S, θ, ϕ)
                Y = vshbasis(Irreducible(), HelicityCovariant(), 0, 0, θ, ϕ, S)
                @test isapproxdefault(Y[0, +1], -√(1/4pi) * WignerD.wignerDjmn(0, 0, 0, 0, θ, ϕ))
                for n in -1:1, basisind in -1:1
                    basisind == 0 && n == 1 && continue
                    @test isapproxdefault(Y[basisind, n], 0)
                end
            end
            for θ in LinRange(0, pi, 10), ϕ in LinRange(0, 2pi, 10)
                cache!(S, θ, ϕ)
                for j in 1:lmax, m in -j:j
                    Y = vshbasis(Irreducible(), HelicityCovariant(), j, m, θ, ϕ, S)
                    @test isapproxdefault(Y[-1, j+1], √(j/8pi) * WignerD.wignerDjmn(j,  1, -m, 0, θ, ϕ))
                    @test isapproxdefault(Y[ 0, j+1], -√((j+1)/4pi) * WignerD.wignerDjmn(j, 0, -m, 0, θ, ϕ))
                    @test isapproxdefault(Y[+1, j+1], √(j/8pi) * WignerD.wignerDjmn(j, -1, -m, 0, θ, ϕ))
                    @test isapproxdefault(Y[-1,  j], √((2j+1)/8pi) * WignerD.wignerDjmn(j,  1, -m, 0, θ, ϕ))
                    @test isapproxdefault(Y[ 0,  j], 0)
                    @test isapproxdefault(Y[+1,  j], -√((2j+1)/8pi) * WignerD.wignerDjmn(j, -1, -m, 0, θ, ϕ))
                    @test isapproxdefault(Y[-1, j-1],  √((j+1)/8pi) * WignerD.wignerDjmn(j,  1, -m, 0, θ, ϕ))
                    @test isapproxdefault(Y[ 0, j-1],  √(j/4pi) * WignerD.wignerDjmn(j, 0, -m, 0, θ, ϕ))
                    @test isapproxdefault(Y[+1, j-1],  √((j+1)/8pi) * WignerD.wignerDjmn(j, -1, -m, 0, θ, ϕ))
                end
            end
        end
    end
    @testset "Polar" begin
        @testset "Hansen" begin
            # l = 0 has only one component
            for θ in LinRange(0, pi, 10), ϕ in LinRange(0, 2pi, 10)
                cache!(S, θ, ϕ, 0)
                Y = vshbasis(Hansen(), Polar(), 0, 0, θ, ϕ, S)
                @test begin
                    res = isapproxdefault(Y[ 1, -1], Ylm[(0,0)])
                    if !res
                        @show θ, ϕ
                    end
                    res
                end
                for n in axes(Y,2), basisind in axes(Y,1)
                    basisind == 1 && n == -1 && continue
                    @test isapproxdefault(Y[basisind, n], 0)
                end
            end
            for θ in LinRange(0, pi, 10), ϕ in LinRange(0, 2pi, 10)
                cache!(S, θ, ϕ)
                for l in 1:lmax, m in -l:l
                    norml⁻ = 1/2*√(((l-m)*(l+m+1))/(l*(l+1)))
                    norml⁺ = 1/2*√(((l+m)*(l-m+1))/(l*(l+1)))
                    Y = vshbasis(Hansen(), Polar(), l, m, θ, ϕ, S)
                    @test isapproxdefault(Y[1, +1],  0)
                    @test isapproxdefault(Y[2, +1],  (m < l ? norml⁻*cis(-ϕ)*Ylm[(l, m+1)] : 0) -
                                          (m > -l ? norml⁺*cis(ϕ)*Ylm[(l, m-1)] : 0)
                        )
                    if sin(θ) != 0
                        @test isapproxdefault(Y[3, +1],  im*m/√(l*(l+1))*1/sin(θ)*Ylm[(l,m)])
                    end
                    @test isapproxdefault(Y[1,  0],  0)
                    if sin(θ) != 0
                        @test isapproxdefault(Y[2,  0],  -m/√(l*(l+1))*1/sin(θ)*Ylm[(l,m)])
                    end
                    @test isapproxdefault(Y[3, 0],  (m < l ? -im*norml⁻*cis(-ϕ)*Ylm[(l, m+1)] : 0) +
                                          (m > -l ? im*norml⁺*cis(ϕ)*Ylm[(l, m-1)] : 0)
                        )
                    @test isapproxdefault(Y[1, -1],  Ylm[(l,m)])
                    @test isapproxdefault(Y[2, -1],  0)
                    @test isapproxdefault(Y[3, -1],  0)

                    if m != 0 && θ == 0
                        # n.Y⁽¹⁾ and n.Y⁽⁰⁾ are zero, and n.Y⁽⁻¹⁾ == YLM which is zero at poles for M != 0
                        @test all(x -> isapproxdefault(x, 0), Y[1, :])
                    end
                end
            end
        end
        @testset "Irreducible" begin
            @testset "Pole" begin
                θ = 0; ϕ = 0
                cache!(S, θ, ϕ)
                for l in 1:lmax, m in -l:l
                    Y = vshbasis(Irreducible(), Polar(), l, m, θ, ϕ, S)
                    # n.Yᴶ⁺¹ and n.Yᴶ⁻¹ ∝ YLM which is zero at poles for M != 0. n.Yᴶ is zero
                    m == 0 && continue
                    @test all(x -> isapproxdefault(x, 0), Y[1, :])
                end
            end
            @testset "J = 0, M = 0" begin
                for θ in LinRange(0, pi, 10), ϕ in LinRange(0, 2pi, 10)
                    cache!(S, θ, ϕ)
                    Y = vshbasis(Irreducible(), Polar(), 0, 0, θ, ϕ, S)
                    @test all(iszero, Y[:,-1:0])
                    @test Y[1, 1] ≈ -1/√(4pi)
                    @test begin
                        res = all(x -> isapproxdefault(0, x), Y[2:3, 1])
                        if !res
                            @show θ, ϕ
                        end
                        res
                    end
                end
            end
        end
    end
    @testset "SphericalCovariant" begin
        @testset "Irreducible" begin
            # For l = 0 only Y₀₀¹ is non-zero
            for θ in LinRange(0, pi, 10), ϕ in LinRange(0, 2pi, 10)
                cache!(S, θ, ϕ, 0)
                Y = vshbasis(Irreducible(), SphericalCovariant(), 0, 0, θ, ϕ, S)
                @test all(iszero, Y[:, -1:0])
                @test isapproxdefault(Y[1, 1], 1/√3 * Ylm[(1,-1)])
                @test isapproxdefault(Y[0, 1], -1/√3 * Ylm[(1, 0)])
                @test isapproxdefault(Y[-1, 1], 1/√3 * Ylm[(1, 1)])
            end
            for θ in LinRange(0, pi, 10), ϕ in LinRange(0, 2pi, 10)
                cache!(S, θ, ϕ)
                for j in 1:lmax, m in -j:j
                    Y = vshbasis(Irreducible(), SphericalCovariant(), j, m, θ, ϕ, S)
                    @test isapproxdefault(Y[ 1, j+1],  √((j-m+1)*(j-m+2)/(2*(j+1)*(2j+3))) * Ylm[(j+1,m-1)])
                    @test isapproxdefault(Y[ 0, j+1], -√((j-m+1)*(j+m+1)/((j+1)*(2j+3)))   * Ylm[(j+1,m)])
                    @test isapproxdefault(Y[-1, j+1],  √((j+m+1)*(j+m+2)/(2*(j+1)*(2j+3))) * Ylm[(j+1,m+1)])
                    @test isapproxdefault(Y[ 1, j], m > -j ? -√((j+m)*(j-m+1)/(2j*(j+1))) * Ylm[(j,m-1)] : 0)
                    @test isapproxdefault(Y[ 0, j],  m/√(j*(j+1)) * Ylm[(j,m)])
                    @test isapproxdefault(Y[-1, j], m < j ? √((j-m)*(j+m+1)/(2j*(j+1))) * Ylm[(j,m+1)] : 0)
                    @test isapproxdefault(Y[ 1, j-1], abs(m - 1) <= j-1 ? √((j+m)*(j+m-1)/(2j*(2j-1))) * Ylm[(j-1,m-1)] : 0)
                    @test isapproxdefault(Y[ 0, j-1], abs(m) <= j-1 ? √((j-m)*(j+m)/(j*(2j-1))) * Ylm[(j-1,m)] : 0)
                    @test isapproxdefault(Y[-1, j-1], abs(m+1) <= j-1 ? √((j-m)*(j-m-1)/(2j*(2j-1))) * Ylm[(j-1,m+1)] : 0)
                end
            end
            @testset "Pole" begin
                @testset "NorthPole" begin
                    cache!(S, NorthPole(), 0, lmax)
                    for j in 1:lmax, m in -j:j
                        # test that θ = 0 and θ = NorthPole() lead to identical results
                        Y0 = vshbasis(Irreducible(), SphericalCovariant(), j, m, 0, 0)
                        Y = vshbasis(Irreducible(), SphericalCovariant(), j, m, NorthPole(), 0, S)
                        @test isapproxdefault(Y0, Y)
                        if abs(m) > 1
                            @test all(iszero, Y)
                        elseif m == 1 || m == -1
                            @test Y[m,j+1] ≈ √(j/8pi)
                            @test Y[m,j] ≈ -m*√((2j+1)/8pi)
                            @test Y[m,j-1] ≈ √((j+1)/8pi)
                        elseif m == 0
                            @test Y[0,j+1] ≈ -√((j+1)/4pi)
                            @test isapproxdefault(Y[0,j], 0)
                            @test Y[0,j-1] ≈ √(j/4pi)
                        end

                        Y = vshbasis(Hansen(), SphericalCovariant(), j, m, NorthPole(), 0, S)
                        if abs(m) > 1
                            @test all(iszero, Y)
                        elseif m == 1 || m == -1
                            @test Y[m, 1] ≈ √((2j+1)/8pi)
                            @test Y[m, 0] ≈ -m*√((2j+1)/8pi)
                            @test isapproxdefault(Y[m,-1], 0)
                        elseif m == 0
                            @test isapproxdefault(Y[0, 1], 0)
                            @test isapproxdefault(Y[0, 0], 0)
                            @test Y[0,-1] ≈ √((2j+1)/4pi)
                        end
                    end
                end
                @testset "SouthPole" begin
                    cache!(S, SouthPole(), 0, lmax)
                    for j in 1:lmax, m in -j:j
                        # test that θ = pi and θ = SouthPole() lead to identical results
                        Y0 = vshbasis(Irreducible(), SphericalCovariant(), j, m, pi, 0)
                        Y = vshbasis(Irreducible(), SphericalCovariant(), j, m, SouthPole(), 0, S)
                        @test isapproxdefault(Y0, Y)
                        if abs(m) > 1
                            @test all(iszero, Y)
                        elseif m == 1 || m == -1
                            @test Y[m,j+1] ≈ (-1)^(j-1)*√(j/8pi)
                            @test Y[m,j] ≈ -(-1)^(j)*m*√((2j+1)/8pi)
                            @test Y[m,j-1] ≈ (-1)^(j-1)*√((j+1)/8pi)
                        elseif m == 0
                            @test Y[0,j+1] ≈ -(-1)^(j-1)*√((j+1)/4pi)
                            @test isapproxdefault(Y[0,j], 0)
                            @test Y[0,j-1] ≈ (-1)^(j-1)*√(j/4pi)
                        end

                        Y = vshbasis(Hansen(), SphericalCovariant(), j, m, SouthPole(), 0, S)
                        if abs(m) > 1
                            @test all(iszero, Y)
                        elseif m == 1 || m == -1
                            @test Y[m, 1] ≈ (-1)^(j-1)*√((2j+1)/8pi)
                            @test Y[m, 0] ≈ -(-1)^(j)*m*√((2j+1)/8pi)
                            @test isapproxdefault(Y[m,-1], 0)
                        elseif m == 0
                            @test isapproxdefault(Y[0, 1], 0)
                            @test isapproxdefault(Y[0, 0], 0)
                            @test Y[0,-1] ≈ (-1)^(j-1)*√((2j+1)/4pi)
                        end
                    end
                end
            end
        end
        @testset "Hansen" begin
            for θ in LinRange(0, pi, 10), ϕ in LinRange(0, 2pi, 10)
                cache!(S, θ, ϕ, 0)
                Y = vshbasis(Hansen(), SphericalCovariant(), 0, 0, θ, ϕ, S)
                @test all(iszero, Y[:, 0:1])
                @test isapproxdefault(Y[ 1, -1], -1/√3 * Ylm[(1,-1)])
                @test isapproxdefault(Y[ 0, -1], 1/√3 * Ylm[(1, 0)])
                @test isapproxdefault(Y[-1, -1], -1/√3 * Ylm[(1, 1)])
            end
            for θ in LinRange(0, pi, 10), ϕ in LinRange(0, 2pi, 10)
                cache!(S, θ, ϕ)
                for l in 1:lmax, m in -l:l
                    Y = vshbasis(Hansen(), SphericalCovariant(), l, m, θ, ϕ, S)
                    @test isapproxdefault(Y[ 1, 1],  √((l-m+1)*(l-m+2)*l/(2*(l+1)*(2l+1)*(2l+3))) * Ylm[(l+1,m-1)] +
                                                    (abs(m - 1) <= l-1 ? √((l+m)*(l+m-1)*(l+1)/(2l*(2l-1)*(2l+1))) * Ylm[(l-1,m-1)] : 0))
                    @test isapproxdefault(Y[ 0, 1], -√((l-m+1)*(l+m+1)*l/((l+1)*(2l+1)*(2l+3)))   * Ylm[(l+1,m)] +
                                                    (abs(m) <= l-1 ? √((l-m)*(l+m)*(l+1)/(l*(2l-1)*(2l+1))) * Ylm[(l-1,m)] : 0))
                    @test isapproxdefault(Y[-1, 1],  √((l+m+1)*(l+m+2)*l/(2*(l+1)*(2l+1)*(2l+3))) * Ylm[(l+1,m+1)] +
                                                    (abs(m+1) <= l-1 ? √((l-m)*(l-m-1)*(l+1)/(2l*(2l-1)*(2l+1))) * Ylm[(l-1,m+1)] : 0))
                    @test isapproxdefault(Y[ 1, 0], m > -l ? -√((l+m)*(l-m+1)/(2l*(l+1))) * Ylm[(l,m-1)] : 0)
                    @test isapproxdefault(Y[ 0, 0],  m/√(l*(l+1)) * Ylm[(l,m)])
                    @test isapproxdefault(Y[-1, 0], m < l ? √((l-m)*(l+m+1)/(2l*(l+1))) * Ylm[(l,m+1)] : 0)
                    @test isapproxdefault(Y[ 1,-1], (abs(m - 1) <= l-1 ? √((l+m)*(l+m-1)/(2*(2l-1)*(2l+1))) * Ylm[(l-1,m-1)] : 0) -
                                                    √((l-m+1)*(l-m+2)/(2*(2l+1)*(2l+3))) * Ylm[(l+1,m-1)])
                    @test isapproxdefault(Y[ 0,-1], (abs(m) <= l-1 ? √((l-m)*(l+m)/((2l-1)*(2l+1))) * Ylm[(l-1,m)] : 0) +
                                                    √((l-m+1)*(l+m+1)/((2l+1)*(2l+3))) * Ylm[(l+1,m)])
                    @test isapproxdefault(Y[-1,-1], (abs(m+1) <= l-1 ? √((l-m)*(l-m-1)/(2*(2l-1)*(2l+1))) * Ylm[(l-1,m+1)] : 0) -
                                                    √((l+m+1)*(l+m+2)/(2*(2l+1)*(2l+3))) * Ylm[(l+1,m+1)])
                end
            end
        end
    end
    @testset "matrix and vectors" begin
        for θ in LinRange(0, pi, 10), ϕ in LinRange(0, 2pi, 10)
            cache!(S, θ, ϕ)
            for YT in [PB(), Irreducible(), Hansen()],
                B in [Cartesian(), Polar(), SphericalCovariant(), HelicityCovariant()]

                for l in 0:lmax, m in -l:l
                    Y = vshbasis(YT, B, l, m, θ, ϕ, S)
                    Yp = parent(Y)
                    for n in axes(Y, 2)
                        Yn = vshbasis(YT, B, l, m, n, θ, ϕ, S)
                        @test begin
                            res = parent(Y[:, n]) ≈ parent(Yn)
                            if !res
                                @show YT, B, l, m, n
                            end
                            res
                        end
                    end
                end
            end
        end

        @testset "genspharm" begin
            for θ in LinRange(0, pi, 10), ϕ in LinRange(0, 2pi, 10)
                cache!(S, θ, ϕ)
                G2 = genspharm(ML(0:lmax), θ, ϕ, S)
                for l in 0:lmax, m in -l:l
                    Y = vshbasis(PB(), HelicityCovariant(), l, m, θ, ϕ, S)
                    G = genspharm(l, m, θ, ϕ, S)
                    @test isapproxdefault(G2[(l,m)], G)
                    nmax = min(1, l)
                    for n in -nmax:nmax
                        Ylmn = genspharm2(l, m, n, θ, ϕ)
                        @test isapproxdefault(Y[n,n], G[n])
                        @test isapproxdefault(G[n], Ylmn)
                    end
                end
            end
        end
    end
end

@testset "one vs multiple" begin
    lmax = 2
    modes = LM(0:lmax)
    modes_high = LM(2:lmax, 2:2)
    S2 = VectorSphericalHarmonics.cache(Float64, lmax)
    for YT in [Irreducible(), PB(), Hansen()], B in [Cartesian(), Polar(), SphericalCovariant(), HelicityCovariant()]
        @testset "$YT $B" begin
            for m in [modes, modes_high], (θ, ϕ) in ((pi/3, pi/3), (0,0))
                cache!(S2, θ, ϕ)
                Y_multiple = vshbasis(YT, B, m, θ, ϕ, S2)
                Y_multiple2 = vshbasis(YT, B, m, θ, ϕ)
                for (ind, (j, m)) in enumerate(m)
                    Y = vshbasis(YT, B, j, m, θ, ϕ, S2)
                    @test isapproxdefault(Y_multiple[ind], Y)
                    @test isapproxdefault(Y_multiple2[ind], Y)
                end
            end
        end
    end

    @testset "genspharm" begin
        for m in [modes, modes_high], (θ, ϕ) in ((pi/3, pi/3), (0,0))
            cache!(S2, θ, ϕ)
            Y_multiple = genspharm(m, θ, ϕ, S2)
            Y_multiple2 = genspharm(m, θ, ϕ)
            for (ind, (j, m)) in enumerate(m)
                Y = genspharm(j, m, θ, ϕ, S2)
                @test isapproxdefault(Y_multiple[ind], Y)
                @test isapproxdefault(Y_multiple2[ind], Y)
            end
        end
    end
end

@testset "Orthogonal" begin
    lmax = 10
    @testset "Locally orthogonal" begin
        @testset "PB" begin
            # The PB basis vectors are directed along the Helicity basis vectors, so they are orthogonal at each point
            for θ in LinRange(0, pi, 10), ϕ in LinRange(0, 2pi, 10)
                cache!(S, θ, ϕ);
                for l in 0:lmax, m in -l:l, B in [SphericalCovariant(), Polar(), HelicityCovariant(), Cartesian()]
                    Y = vshbasis(PB(), B, l, m, θ, ϕ, S)
                    Yp = parent(Y)
                    M = Yp'*Yp
                    for j in axes(M, 2), i in axes(M, 1)
                        i == j && continue
                        @test isapprox(M[i, j], 0, atol=1e-13, rtol=1e-8)
                    end
                end
            end
        end
        @testset "Hansen" begin
            # For Hansen basis, the local orthogonality uses transpose instead of adjoint.
            for θ in LinRange(0, pi, 10), ϕ in LinRange(0, 2pi, 10)
                cache!(S, θ, ϕ);
                for l in 0:lmax, m in -l:l
                    for B in [Polar(), Cartesian()]
                        Y = vshbasis(Hansen(), B, l, m, θ, ϕ, S)
                        Yp = parent(Y)
                        M = transpose(Yp)*Yp
                        for j in axes(M, 2), i in axes(M, 1)
                            i == j && continue
                            @test begin
                                res = isapprox(M[i, j], 0, atol=1e-13, rtol=1e-8)
                                if !res
                                    @show B, l, m, θ, ϕ
                                end
                                res
                            end
                        end
                    end
                    # For the complex bases, we need to use eμ⋅eν = (-1)^μ δμ,-ν
                    for B in [SphericalCovariant(), HelicityCovariant()]
                        Y = vshbasis(Hansen(), B, l, m, θ, ϕ, S)
                        for j in axes(Y, 2), i in axes(Y, 2)
                            i == j && continue
                            M = sum(Y[b,j] * Y[-b,i] * (-1)^b for b = -1:1)
                            @test begin
                                res = isapprox(M, 0, atol=1e-13, rtol=1e-8)
                                if !res
                                    @show B, l, m, θ, ϕ
                                end
                                res
                            end
                        end
                    end
                end
            end
        end
    end
    @testset "Orthonormal" begin
        lmax = 3
        function f(θ, ϕ, YT, B, l, m)
            cache!(S, θ, ϕ)
            Y = vshbasis(YT, B, l, m, θ, ϕ, S)
            Yp = parent(Y)
            (Yp'*Yp) .* sin(θ)
        end
        @testset "PB" begin
            for l in 0:lmax, m in -l:l, B in [SphericalCovariant(), Polar(), HelicityCovariant(), Cartesian()]
                M = hcubature(x -> f(x..., PB(), B, l, m), [0, 0], [pi, 2pi])[1]
                if l == 0
                    for I in CartesianIndices(M)
                        i, j = Tuple(I)
                        i == 2 && j == 2 && continue
                        @test isapproxdefault(M[i,j], 0)
                    end
                    @test M[2,2] ≈ 1
                else
                    @test M ≈ I
                end
            end
        end
        @testset "Irreducible" begin
            for l in 0:lmax, m in -l:l, B in [SphericalCovariant(), Polar(), HelicityCovariant(), Cartesian()]
                M = hcubature(x -> f(x..., Irreducible(), B, l, m), [0, 0], [pi, 2pi])[1]
                if l == 0
                    for I in CartesianIndices(M)
                        i, j = Tuple(I)
                        i == 3 && j == 3 && continue
                        @test isapproxdefault(M[i,j], 0)
                    end
                    @test M[3,3] ≈ 1
                else
                    @test M ≈ I
                end
            end
        end
        @testset "Hansen" begin
            for l in 0:lmax, m in -l:l, B in [SphericalCovariant(), Polar(), HelicityCovariant(), Cartesian()]
                M = hcubature(x -> f(x..., Hansen(), B, l, m), [0, 0], [pi, 2pi])[1]
                if l == 0
                    for I in CartesianIndices(M)
                        i, j = Tuple(I)
                        i == 1 && j == 1 && continue
                        @test isapproxdefault(M[i,j], 0)
                    end
                    @test M[1,1] ≈ 1
                else
                    @test M ≈ I
                end
            end
        end
    end
end

@testset "Addition theorem" begin
    lmax = 10
    @testset "same point" begin
        @testset "Irreducible" begin
            θ, ϕ = pi/3, pi/3
            for B in [Polar(), SphericalCovariant(), Cartesian(), HelicityCovariant()]
                Y1_all = vshbasis(Irreducible(), B, ML(0:lmax), θ, ϕ);
                for j in 1:lmax, l1 in j-1:j+1, l2 in j-1:j+1
                    S = zero(Float64)
                    for m in -j:j
                        Yjm = Y1_all[(j,m)]
                        Yjml1 = parent(Yjm[:,l1])
                        Yjml2 = parent(Yjm[:,l2])
                        S += dot(Yjml1, Yjml2)
                    end
                    if l1 == l2
                        @test isapproxdefault(S, (2j+1)/4pi)
                    else
                        @test isapproxdefault(S, 0)
                    end
                end
            end
        end
        @testset "Hansen" begin
            θ, ϕ = pi/3, pi/3
            for B in [Polar(), SphericalCovariant(), Cartesian(), HelicityCovariant()]
                Y1_all = vshbasis(Hansen(), B, ML(0:lmax), θ, ϕ);
                for j in 1:lmax, λ1 in -1:+1, λ2 in -1:+1
                    S = zero(Float64)
                    for m in -j:j
                        Yjm = Y1_all[(j,m)]
                        Yjmλ1 = parent(Yjm[:,λ1])
                        Yjmλ2 = parent(Yjm[:,λ2])
                        S += dot(Yjmλ1, Yjmλ2)
                    end
                    if λ1 == λ2
                        @test isapproxdefault(S, (2j+1)/4pi)
                    else
                        @test isapproxdefault(S, 0)
                    end
                end
            end
        end
    end
    @testset "different points" begin
        @testset "Irreducible" begin
            θ1, ϕ1 = pi/3, pi/3
            θ2, ϕ2 = pi/4, pi/4
            cosω12 = cos(θ1)cos(θ2) + sin(θ1)sin(θ2)cos(ϕ1 - ϕ2)
            P = LegendrePolynomials.collectPl(cosω12, lmax = lmax + 1);
            allmodes = ML(0:lmax)
            # We only test in bases that are constant, so they are identical at both the points
            for B in [SphericalCovariant(), Cartesian()]
                Y1_all = vshbasis(Irreducible(), B, allmodes, θ1, ϕ1);
                Y2_all = vshbasis(Irreducible(), B, allmodes, θ2, ϕ2);
                for j in 1:lmax, l1 in j-1:j+1, l2 in j-1:j+1
                    x = zero(Float64)
                    for m in -j:j
                        Y1jm = Y1_all[(j,m)]
                        Y2jm = Y2_all[(j,m)]
                        Y1jml1 = parent(Y1jm[:,l1])
                        Y2jml2 = parent(Y2jm[:,l2])
                        x += dot(Y1jml1, Y2jml2)
                    end
                    @test isapproxdefault(imag(x), 0)
                    if l1 == l2
                        @test begin
                            res = isapproxdefault(x, (2j+1)/4pi * P[l1])
                            if !res
                                @show j, l1, B, x, (2j+1)/4pi * P[l1]
                            end
                            res
                        end
                    else
                        @test begin
                            res = isapproxdefault(x, 0)
                            if !res
                                @show j, l1, l2, B, x
                            end
                            res
                        end
                    end
                end
            end
        end
    end
end

@testset "Complex conjugation" begin
    @testset "Irreducible" begin
        # For l = 0 Y₀₀¹ is purely real
        for θ in LinRange(0, pi, 10), ϕ in LinRange(0, 2pi, 10)
            cache!(S, θ, ϕ, 0)
            Y = vshbasis(Irreducible(), Polar(), 0, 0, θ, ϕ, S)
            @test all(x -> isapproxdefault(imag(x), 0), Y)
        end
        for θ in LinRange(0, pi, 10), ϕ in LinRange(0, 2pi, 10)
            cache!(S, θ, ϕ)
            for j in 1:5, m in -j:j
                Ym = vshbasis(Irreducible(), Polar(), j, m, θ, ϕ, S)
                Ymconj = conj.(Ym)
                Y₋m = vshbasis(Irreducible(), Polar(), j, -m, θ, ϕ, S)
                for l in axes(Ym, 2)
                    Ymnconj = Ymconj[:, l]
                    Y₋mn = Y₋m[:, l]
                    phase = (-1)^(j + m + l + 1)
                    @test all(isapproxdefault, zip(Ymnconj, phase .* Y₋mn))
                end
            end
        end
    end
    @testset "Hansen" begin
        # For l = 0 Y₀₀⁽⁻¹⁾ is purely real
        for θ in LinRange(0, pi, 10), ϕ in LinRange(0, 2pi, 10)
            cache!(S, θ, ϕ, 0)
            Y = vshbasis(Hansen(), Polar(), 0, 0, θ, ϕ, S)
            @test all(x -> isapproxdefault(imag(x), 0), Y)
        end
        for θ in LinRange(0, pi, 10), ϕ in LinRange(0, 2pi, 10)
            cache!(S, θ, ϕ)
            for l in 1:5, m in -l:l
                Ym = vshbasis(Hansen(), Polar(), l, m, θ, ϕ, S)
                Ymconj = conj.(Ym)
                Y₋m = vshbasis(Hansen(), Polar(), l, -m, θ, ϕ, S)
                for λ in axes(Ym, 2)
                    Ymnconj = Ymconj[:, λ]
                    Y₋mn = Y₋m[:, λ]
                    phase = (-1)^(m + λ + 1)
                    @test all(isapproxdefault.(Ymnconj, phase .* Y₋mn))
                end
            end
        end
    end
    @testset "PB" begin
        for θ in LinRange(0, pi, 10), ϕ in LinRange(0, 2pi, 10)
            cache!(S, θ, ϕ)
            for j in 1:5, m in -j:j
                Ym = vshbasis(PB(), Polar(), j, m, θ, ϕ, S)
                Ymconj = conj.(Ym)
                Y₋m = vshbasis(PB(), Polar(), j, -m, θ, ϕ, S)
                for n in axes(Ym, 2)
                    Y₋m₋n = Y₋m[:, -n]
                    Ymn = Ym[:, n]
                    phase = (-1)^m
                    @test begin
                        res = all(isapproxdefault, zip(conj.(Y₋m₋n), phase .* Ymn))
                        if !res
                            @show j, m, n, θ, ϕ
                        end
                        res
                    end
                end
            end
        end
        for θ in LinRange(0, pi, 10), ϕ in LinRange(0, 2pi, 10)
            cache!(S, θ, ϕ)
            for j in 1:5, m in -j:j
                Ym = vshbasis(PB(), HelicityCovariant(), j, m, θ, ϕ, S)
                Ymconj = conj.(Ym)
                Y₋m = vshbasis(PB(), HelicityCovariant(), j, -m, θ, ϕ, S)
                for n in axes(Ym, 2)
                    Y₋m₋n = Y₋m[-n, -n]
                    Ymn = Ym[n, n]
                    phase = (-1)^(m + n)
                    @test begin
                        res = all(isapproxdefault, zip(conj.(Y₋m₋n), phase .* Ymn))
                        if !res
                            @show j, m, n, θ, ϕ
                        end
                        res
                    end
                end
            end
        end
    end
end

@testset "Coordinate inversion" begin
    jmax = 4
    @testset "Irreducible" begin
        for θ in LinRange(0, pi, 10), ϕ in LinRange(0, 2pi, 10)
            cache!(S, θ, ϕ, jmax)
            for j in 0:jmax, m in -j:j
                for B in [Cartesian(), SphericalCovariant()]
                    Yjm = vshbasis(Irreducible(), B, j, m, θ, ϕ, S)
                    Yjmrot = vshbasis(Irreducible(), B, j, m, pi - θ, pi + ϕ)
                    for l in axes(Yjm, 2)
                        Yjml = Yjm[:, l] * (-1)^(l+1)
                        Yjmlrot = -Yjmrot[:, l]
                        @test begin
                            res = all(isapproxdefault, zip(Yjml, Yjmlrot))
                            if !res
                                @show j, m, θ, ϕ
                            end
                            res
                        end
                    end
                end
            end
        end
    end
    @testset "Hansen" begin
        for θ in LinRange(0, pi, 10), ϕ in LinRange(0, 2pi, 10)
            cache!(S, θ, ϕ, jmax)
            for j in 0:4, m in -j:j
                for B in [Cartesian(), SphericalCovariant()]
                    Yjm = vshbasis(Hansen(), B, j, m, θ, ϕ, S)
                    Yjmrot = vshbasis(Hansen(), B, j, m, pi - θ, pi + ϕ)
                    for λ in axes(Yjm, 2)
                        Yjml = Yjm[:, λ] * (-1)^(j + λ +1)
                        Yjmlrot = -Yjmrot[:, λ]
                        @test all(isapproxdefault, zip(Yjml, Yjmlrot))
                    end
                end
            end
        end
    end
end

@testset "Coordinate rotation" begin
    lmax = 3
    modes = ML(0:lmax)

    S2 = VectorSphericalHarmonics.cache(0, 0, lmax);
    S′ = VectorSphericalHarmonics.cache(0, 0, lmax);
    Dvec = OffsetArray([zeros(ComplexF64, 2l+1, 2l+1) for l in 0:lmax], 0:lmax);

    @testset "arbitrary point" begin
        cartvec(θ, ϕ) = SVector{3}(sin(θ)cos(ϕ), sin(θ)sin(ϕ), cos(θ))
        function polcoords(n)
            nx, ny, nz = normalize(n)
            θ = acos(nz)
            if abs(nz) == 1
                # ϕ is not defined uniquely at the poles
                return promote(θ, zero(θ))
            else
                ϕ = atan(ny, nx)
                return promote(θ, ϕ)
            end
        end

        @testset "genspharm" begin
            for θ in LinRange(0, pi, 6), ϕ in LinRange(0, 2pi, 6)
                VectorSphericalHarmonics.cache!(S2, θ, ϕ);
                Y = genspharm(modes, θ, ϕ, S2);
                n = cartvec(θ, ϕ)
                Un = basisconversionmatrix(Cartesian(), HelicityCovariant(), θ, ϕ);

                for θ′ in LinRange(0, pi, 6), ϕ′ in LinRange(0, 2pi, 6)
                    VectorSphericalHarmonics.cache!(S′, θ′, ϕ′);
                    α, β, γ = ϕ, θ-θ′, -ϕ′
                    R = RotZYZ(α, β, γ)
                    R⁻¹ = inv(R)
                    # this specific rotation satisfies the condition U′n * R' * Un' == I,
                    # so the basis rotation matrix may be left out
                    # This happens because Un for polar == RotYZ(-θ, -ϕ) (albeit with permuted rows), so
                    # U′n * R' * Un' = RotYZ(-θ′, -ϕ′) * RotZYZ(ϕ′, θ′ - θ, -ϕ) * RotZY(ϕ, θ)  ≈ I
                    # Un for HelicityCovariant is related to that for the polar matrix through a unitary
                    # transformation, so the same relation holds
                    U′n = basisconversionmatrix(Cartesian(), HelicityCovariant(), θ′, ϕ′);
                    @test U′n * R⁻¹ * Un' ≈ I
                    for l in 0:lmax
                        Dp = wignerD!(Dvec[l], l, α, β, γ);
                        D = OffsetArray(Dp, -l:l, -l:l)
                        for m in -l:l
                            Ylθϕrot_m = sum(D[m′, m] * Y[(l, m′)] for m′ in -l:l)
                            Ylmθ′ϕ′ = genspharm(l, m, θ′, ϕ′, S′)
                            @test isapprox(Ylmθ′ϕ′, Ylθϕrot_m, atol = 1e-13, rtol = 1e-8)
                        end
                    end
                end

                # obtain the second point from an arbitrary rotation
                # In general we need to retain the basis conversion matrix to compute the components correctly
                for α in LinRange(0, 2pi, 6), β in LinRange(0, pi, 6), γ in LinRange(0, 2pi, 6)
                    # the rotation that transforms between frames S′ = RS
                    R = RotZYZ(α, β, γ)
                    R⁻¹ = inv(R)
                    n′ = R⁻¹ * n
                    θ′, ϕ′ = polcoords(n′)
                    VectorSphericalHarmonics.cache!(S′, θ′, ϕ′);
                    U′n = basisconversionmatrix(Cartesian(), HelicityCovariant(), θ′, ϕ′);
                    M⁻¹ = U′n * R⁻¹ * Un'
                    @test M⁻¹ ≈ Diagonal(M⁻¹)
                    for l in 0:lmax
                        Dp = wignerD!(Dvec[l], l, α, β, γ);
                        D = OffsetArray(Dp, -l:l, -l:l)
                        for m in -l:l
                            Ylθϕrot_m = M⁻¹ * sum(D[m′, m] * parent(Y[(l, m′)]) for m′ in -l:l)
                            Ylmθ′ϕ′ = parent(genspharm(l, m, θ′, ϕ′, S′))
                            @test isapprox(Ylmθ′ϕ′, Ylθϕrot_m, atol = 1e-13, rtol = 1e-8)
                        end
                    end
                end
            end
        end

        @testset "VSH" begin
            for θ in LinRange(0, pi, 6), ϕ in LinRange(0, 2pi, 6)
                n = cartvec(θ, ϕ)
                VectorSphericalHarmonics.cache!(S2, θ, ϕ);

                for YT in [Irreducible(), Hansen(), PB()], B in [Cartesian(), SphericalCovariant(), HelicityCovariant(), Polar()]
                    Un = basisconversionmatrix(Cartesian(), B, θ, ϕ)
                    Y = vshbasis(YT, B, modes, θ, ϕ, S2);

                    # We evaluate the rotation of the dot product between VSH and a vector field
                    # In this case we choose a constant vector field directed along x̂ at each point,
                    # where x̂ is directed along the x axis of the frame S (which differs from that of S′ in general)
                    # We evaluate the components in the basis B at both points
                    x_n = Un * SVector{3}(1,0,0)

                    if B ∈ (Polar(), HelicityCovariant())
                        for θ′ in LinRange(0, pi, 6), ϕ′ in LinRange(0, 2pi, 6)
                            VectorSphericalHarmonics.cache!(S′, θ′, ϕ′);
                            α, β, γ = ϕ, θ-θ′, -ϕ′
                            R = RotZYZ(α, β, γ)
                            # this specific rotation satisfies the condition U′n * R' * Un' == I for the
                            # Polar and HelicityCovariant bases, so the basis rotation matrix may be left out
                            # This happens because Un for polar == RotYZ(-θ, -ϕ) (albeit with permuted rows), so
                            # U′n * R' * Un' = RotYZ(-θ′, -ϕ′) * RotZYZ(ϕ′, θ′ - θ, -ϕ) * RotZY(ϕ, θ) ≈ I
                            # Un for HelicityCovariant is related to that for the polar matrix through a unitary
                            # transformation, so the same relation holds
                            U′n = basisconversionmatrix(Cartesian(), B, θ′, ϕ′);
                            if B === Polar()
                                # We go from (r,θ,ϕ) to (θ,ϕ,r)
                                @test Un[[2,3,1],:] ≈ RotYZ(-θ, -ϕ)
                                @test U′n[[2,3,1],:] ≈ RotYZ(-θ′, -ϕ′)
                            end
                            @test U′n * R' * Un' ≈ I
                            for l in 0:lmax
                                Dp = wignerD!(Dvec[l], l, α, β, γ);
                                D = OffsetArray(Dp, -l:l, -l:l)
                                for m in -l:l
                                    Ylθϕrot_m = sum(D[m′, m] * Y[(l, m′)] for m′ in -l:l)
                                    Ylmθ′ϕ′ = vshbasis(YT, B, l, m, θ′, ϕ′, S′)
                                    @test isapprox(Ylmθ′ϕ′, Ylθϕrot_m, atol = 1e-13, rtol = 1e-8)
                                end
                            end
                        end
                    end

                    # obtain the second point from an arbitrary rotation
                    # In general we need to retain the basis conversion matrix to compute the components correctly
                    for α in LinRange(0, 2pi, 6), β in LinRange(0, pi, 6), γ in LinRange(0, 2pi, 6)
                        # the rotation that transforms between frames S′ = RS
                        R = RotZYZ(α, β, γ)
                        R⁻¹ = inv(R)
                        n′ = R⁻¹ * n
                        θ′, ϕ′ = polcoords(n′)
                        VectorSphericalHarmonics.cache!(S′, θ′, ϕ′);
                        U′n = basisconversionmatrix(Cartesian(), B, θ′, ϕ′);

                        M⁻¹ = (U′n * R⁻¹ * Un')
                        if B === HelicityCovariant()
                            @test M⁻¹ ≈ Diagonal(M⁻¹)
                        end
                        M = inv(M⁻¹)

                        x_n′ = M⁻¹ * x_n

                        for l in 0:lmax
                            Dp = wignerD!(Dvec[l], l, α, β, γ);
                            D = OffsetArray(Dp, -l:l, -l:l)
                            for m in -l:l
                                Ylθϕrot_m = M⁻¹ * sum(D[m′, m] * parent(Y[(l, m′)]) for m′ in -l:l)
                                Ylmθ′ϕ′ = parent(vshbasis(YT, B, l, m, θ′, ϕ′, S′))
                                @test isapprox(Ylmθ′ϕ′, Ylθϕrot_m, atol = 1e-13, rtol = 1e-8)
                                # (M⁻¹ x)⋅Ylm′(n′) = (M⁻¹ x)⋅M⁻¹∑Dlmm′ Ylm(n) = x⋅∑Dlmm′ Ylm(n)
                                Ylθϕrot_m = x_n' * sum(D[m′, m] * parent(Y[(l, m′)]) for m′ in -l:l)
                                Ylmθ′ϕ′ = x_n′' * parent(vshbasis(YT, B, l, m, θ′, ϕ′, S′))
                                @test isapprox(Ylmθ′ϕ′, Ylθϕrot_m, atol = 1e-13, rtol = 1e-8)
                                # x.Ylm′(n′) = x.M⁻¹∑Dlmm′ Ylm(n) = (Mx)⋅∑Dlmm′ Ylm(n)
                                Ylθϕrot_m = (M*x_n)' * sum(D[m′, m] * parent(Y[(l, m′)]) for m′ in -l:l)
                                Ylmθ′ϕ′ = x_n' * parent(vshbasis(YT, B, l, m, θ′, ϕ′, S′))
                                @test isapprox(Ylmθ′ϕ′, Ylθϕrot_m, atol = 1e-13, rtol = 1e-8)
                            end
                        end
                    end
                end
            end
        end
    end
end

@testset "VSHcache" begin
    modes = LM(0:1, 0:1)
    lmax = 2
    θ1, ϕ1 = pi/3, pi/2
    θ2, ϕ2 = pi/4, pi/6
    @testset "vshbasis!" begin
        V = @inferred VectorSphericalHarmonics.VSHCache(Float64, Irreducible(), Polar(), θ1, ϕ1, modes);
        M = vshbasis(Irreducible(), Polar(), modes, θ1, ϕ1)
        @test getY(V) == M
        @test eltypeY(V) == eltype(getY(V))

        for θ in Any[θ2, NorthPole(), SouthPole()]
            vshbasis!(V, Irreducible(), Polar(), θ, ϕ2)
            M = vshbasis(Irreducible(), Polar(), modes, θ, ϕ2)
            @test getY(V) == M
        end
    end
    @testset "genspharm!" begin
        V = @inferred VectorSphericalHarmonics.VSHCache(Float64, θ1, ϕ1, modes);
        M = genspharm(modes, θ1, ϕ1)
        @test getY(V) == M

        for θ in Any[θ2, NorthPole(), SouthPole()]
            genspharm!(V, θ, ϕ2)
            M = genspharm(modes, θ, ϕ2)
            @test getY(V) == M
        end
    end
    @testset "uninitialized" begin
        V = VectorSphericalHarmonics.VSHCache(PB(), Cartesian(), LM(0:1))
        @test eltype(VectorSphericalHarmonics.eltypeY(V)) == ComplexF64
        @test all(iszero, VectorSphericalHarmonics.getY(V))

        V = VectorSphericalHarmonics.VSHCache(LM(0:1))
        @test eltype(VectorSphericalHarmonics.eltypeY(V)) == ComplexF64
        @test all(iszero, VectorSphericalHarmonics.getY(V))
    end
end

@testset "broadcast tags" begin
    function components_of_x̂(Basis, (θ,ϕ))
        B = VectorSphericalHarmonics.basisconversionmatrix(Cartesian(), Basis, θ, ϕ)
        B * SVector{3}(1,0,0)
    end
    for B in [Cartesian(), Polar(), SphericalCovariant(), HelicityCovariant()]
        x1, x2 = components_of_x̂.(B, ((pi/2, 0), (pi/2, pi/4)))
        @test x1 == components_of_x̂(B, (pi/2, 0))
        @test x2 == components_of_x̂(B, (pi/2, pi/4))
    end
    V = vshbasis.(Hansen(), Polar(), 1:2, -1, pi/2, 0)
    @test V[1] == vshbasis(Hansen(), Polar(), 1, -1, pi/2, 0)
    @test V[2] == vshbasis(Hansen(), Polar(), 2, -1, pi/2, 0)
    V = vshbasis.(Hansen(), Polar(), 1, -1, (pi/2, pi/4), 0)
    @test V[1] == vshbasis(Hansen(), Polar(), 1, -1, pi/2, 0)
    @test V[2] == vshbasis(Hansen(), Polar(), 1, -1, pi/4, 0)
end

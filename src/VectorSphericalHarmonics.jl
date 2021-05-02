module VectorSphericalHarmonics

using OffsetArrays
using StaticArrays
using LinearAlgebra
using SphericalHarmonics
using SphericalHarmonicModes
using SphericalHarmonicArrays
import SphericalHarmonics: Pole, NorthPole, SouthPole

export vshbasis
export genspharm
export PB
export Irreducible
export Hansen
export SphericalCovariant
export Polar
export HelicityCovariant
export Cartesian

"""
    AbstractVSH

Abstract supertype of vector spherical harmonics
"""
abstract type AbstractVSH end
@doc raw"""
    Irreducible <: AbstractVSH

Vector spherical harmonics that are eigenfunctions of irreducible representationa of total, orbital as well as spin angular momenta.
They may be constructed by coupling scalar spherical harmonics ``Y_{L m}\left(\hat{n}\right)`` with the spherical basis
vectors ``\chi_\mu`` as

```math
\mathbf{Y}_{J M}^L\left(\hat{n}\right) = \sum_{m \mu} C^{J M}_{L m 1 \mu} Y_{L m}\left(\hat{n}\right) \chi_\mu,
```
"""
struct Irreducible <: AbstractVSH end

@doc raw"""
    Hansen <: AbstractVSH

Hansen vector spherical harmonics ``\mathbf{H}_{J M}^{(\lambda)}\left(\hat{n}\right)`` that are related to the [`Irreducible`](@ref)
harmonics through

```math
\begin{aligned}
\mathbf{H}_{JM}^{\left(-1\right)}\left(\hat{n}\right) & =-\sqrt{\frac{J+1}{2J+1}}\mathbf{Y}_{JM}^{J+1}\left(\hat{n}\right)+\sqrt{\frac{J}{2J+1}}\mathbf{Y}_{JM}^{J-1}\left(\hat{n}\right),\\
\mathbf{H}_{JM}^{\left(0\right)}\left(\hat{n}\right) & =\mathbf{Y}_{JM}^{J}\left(\hat{n}\right),\\
\mathbf{H}_{JM}^{\left(1\right)}\left(\hat{n}\right) & =\sqrt{\frac{J}{2J+1}}\mathbf{Y}_{JM}^{J+1}\left(\hat{n}\right)+\sqrt{\frac{J+1}{2J+1}}\mathbf{Y}_{JM}^{J-1}\left(\hat{n}\right),
\end{aligned}
```
"""
struct Hansen <:AbstractVSH end

@doc raw"""
    PB <: AbstractVSH

Phinney-Burridge vector spherical harmonics ``\mathbf{P}_{J M}^\gamma\left(\hat{n}\right)`` that are related to the [`Hansen`](@ref)
harmonics through

```math
\begin{aligned}
\mathbf{P}_{JM}^{+1}\left(\hat{n}\right) & =\frac{1}{\sqrt{2}}\left(\mathbf{H}_{JM}^{\left(1\right)}\left(\hat{n}\right)-\mathbf{H}_{JM}^{\left(0\right)}\left(\hat{n}\right)\right),\\
\mathbf{P}_{JM}^{0}\left(\hat{n}\right) & =\mathbf{H}_{JM}^{\left(-1\right)}\left(\hat{n}\right),\\
\mathbf{P}_{JM}^{-1}\left(\hat{n}\right) & =\frac{1}{\sqrt{2}}\left(\mathbf{H}_{JM}^{\left(1\right)}\left(\hat{n}\right)+\mathbf{H}_{JM}^{\left(0\right)}\left(\hat{n}\right)\right).
\end{aligned}
```
"""
struct PB <: AbstractVSH end

include("basis.jl")

struct VSHCache{T}
    S :: T
end

@doc raw"""
    cache([T::Type = Float64], θ, ϕ, jmax)

Pre-allocate a set of scalar spherical harmonics ``Y_{j m}(\theta,\phi)`` that go into evaluating the vector harmonics.
The cutoff `jmax` corresponds to the maximum degree of vector harmonics that we seek to evaluate.

The type `T` sets the precision used to evaluate the harmonics.
The spherical harmonics evaluated will be of type `Complex{T}`.
"""
cache(θ, ϕ, jmax) = cache(Float64, θ, ϕ, jmax)
function cache(T::Type, θ, ϕ, jmax)
    S = SphericalHarmonics.cache(T, jmax + 1)
    shcache!(S, θ, ϕ)
    return VSHCache(S)
end

@doc raw"""
    cache!(S, θ, ϕ, [jmax])

Update a pre-allocated set of scalar spherical harmonics ``Y_{j m}(\theta,\phi)``.
The cutoff `jmax` corresponds to the maximum degree of vector harmonics that we seek to evaluate.
If this is not provided, the cutoff angular degree in `S` will remain unchanged.
"""
cache!(S::VSHCache, θ, ϕ, jmax...) = (shcache!(S.S, θ, ϕ, jmax...); return S)
function shcache!(S, θ, ϕ, jmax)
    computePlmcostheta!(S, θ, jmax + 1)
    computeYlm!(S, θ, ϕ, jmax + 1)
    return S
end
function shcache!(S, θ, ϕ)
    computePlmcostheta!(S, θ)
    computeYlm!(S, θ, ϕ)
    return S
end

_getY(S::VSHCache) = S.S.Y
_getP(S::VSHCache) = S.S.P
_eltypeY(S::VSHCache) = eltype(_getY(S))
_eltypeP(S::VSHCache) = eltype(_getP(S))

_maybewrapoffset(v, ::Union{Cartesian, Polar}) = v
_maybewrapoffset(v, B::Basis) = OffsetArray(v, _basisinds(B))

_vectorinds(::Irreducible, j) = j .+ (-1:1)
_vectorinds(::Hansen, j) = -1:1
_vectorinds(::PB, j) = -1:1

"""
    vshbasis(Y::AbstractVSH, B::Basis, j::Integer, m::Integer, n::Integer, θ, ϕ, [S = VectorSphericalHarmonics.cache(θ, ϕ, j)])

Evaluate the components of the vector spherical harmonics ``\\mathbf{Y}_{j m}^n(θ, ϕ)`` in the basis `B`.
A pre-allocated array of scalar spherical harmonics `S` may be passed as the final argument.
"""
function vshbasis(YT::AbstractVSH, B::Basis, j, m, n, θ, ϕ, S::VSHCache = cache(θ, ϕ, j))
    Y = vshbasis(YT, B, j, m, θ, ϕ, S)
    Yp = parent(Y)
    ind2 = first(axes(Y,2))
    vs = Yp[:, n - ind2 + 1]
    _maybewrapoffset(vs, B)
end

_PBHelicitycheck(Y, B, M) = M
_PBHelicitycheck(::PB, ::HelicityCovariant, M) = Diagonal(M[SVector{3}(diagind(M))])

"""
    vshbasis(Y::AbstractVSH, B::Basis, j::Integer, m::Integer, θ, ϕ, [S = VectorSphericalHarmonics.cache(θ, ϕ, j)])

Evaluate a set of vector spherical harmonics ``\\mathbf{Y}_{j m}^\\alpha(θ, ϕ)`` for valid values of ``\\alpha``,
and return their components in the basis `B`.
A pre-allocated array of scalar spherical harmonics `S` may be passed as the final argument.
"""
function vshbasis(Y::AbstractVSH, B::Basis, j::Integer, m::Integer, θ, ϕ, S::VSHCache = cache(θ, ϕ, j))
    M = _vshbasis_angle(Y, B, j, m, θ, ϕ, S)
    M2 = _PBHelicitycheck(Y, B, M)
    OffsetArray(M2, _basisinds(B), _vectorinds(Y, j))
end

function _vshbasis_angle(Y::AbstractVSH, B::Basis, j, m, θ, ϕ, S::VSHCache)
    if θ == 0
        return _vshbasis_angle(Y, B, j, m, NorthPole(), ϕ, S)
    elseif θ == pi
        return _vshbasis_angle(Y, B, j, m, SouthPole(), ϕ, S)
    else
        return _vshbasis(Y, B, j, m, θ, ϕ, S)
    end
end

# VSH at poles are zero except for M = 0, ±1
function _vshbasis_angle(Y::AbstractVSH, B::Basis, j, m, θ::Pole, ϕ, S::VSHCache)
    if abs(m) > 1
        T = _eltypeY(S)
        M = SMatrix{3,3,T}(ntuple(x -> zero(T), Val(9)))
    else
        M = _vshbasis(Y, B, j, m, θ, ϕ, S)
    end
    return M
end

function _vshbasis(::PB, B, j, m, θ, ϕ, S)
    H = vshbasis(Hansen(), B, j, m, θ, ϕ, S)
    Hp_static = SMatrix{3,3}(parent(H))
    Hm1 = Hp_static[:, 1]
    H0  = Hp_static[:, 2]
    Hp1 = Hp_static[:, 3]

    T = _eltypeP(S)
    Ym1 = 1/√T(2) * (Hp1 + H0)
    Y0 = Hm1
    Yp1 = 1/√T(2) * (Hp1 - H0)

    M = SMatrix{3,3}((Ym1..., Y0..., Yp1...))
end

function __vshbasis(::Hansen, B, j, m, θ, ϕ, S)
    Y = vshbasis(Irreducible(), B, j, m, θ, ϕ, S)
    Yp_static = SMatrix{3,3}(parent(Y))
    Yjm1 = Yp_static[:, 1]
    Yj  = Yp_static[:, 2]
    Yjp1 = Yp_static[:, 3]

    T = _eltypeP(S)
    Hp1 = √(T(j)/(2j+1)) * Yjp1 + √(T(j+1)/(2j+1)) * Yjm1
    H0 = Yj
    Hm1 = -√(T(j+1)/(2j+1)) * Yjp1 + √(T(j)/(2j+1)) * Yjm1

    M = SMatrix{3,3}((Hm1..., H0..., Hp1...))
    OffsetArray(M, _basisinds(B), _vectorinds(Hansen(), j))
end
function _vshbasis(::Hansen, B, j, m, θ, ϕ, S)
    M = __vshbasis(Hansen(), B, j, m, θ, ϕ, S)
    parent(M)
end
function _vshbasis(::Hansen, B::Polar, j, m, θ, ϕ, S::VSHCache)
    HO = __vshbasis(Hansen(), B, j, m, θ, ϕ, S)
    Hm1 = SVector{3}(HO[1,-1], zero(eltype(HO)), zero(eltype(HO)))
    H0 = SVector{3}(zero(eltype(HO)), HO[2,0], HO[3,0])
    Hp1 = SVector{3}(zero(eltype(HO)), HO[2,1], HO[3,1])
    SMatrix{3,3}(Hm1..., H0..., Hp1...)
end
function _vshbasis(::Hansen, B::HelicityCovariant, j, m, θ, ϕ, S::VSHCache)
    HO = __vshbasis(Hansen(), B, j, m, θ, ϕ, S)
    Hm1 = SVector{3}(zero(eltype(HO)), HO[0,-1], zero(eltype(HO)))
    H0 = SVector{3}(HO[-1,0], zero(eltype(HO)), HO[1,0])
    Hp1 = SVector{3}(HO[-1,1], zero(eltype(HO)), HO[1,1])
    SMatrix{3,3}(Hm1..., H0..., Hp1...)
end

function _vshbasis(::Irreducible, B::SphericalCovariant, j, m, θ, ϕ, S::VSHCache)
    Y = _getY(S)
    TC = eltype(Y)
    T = _eltypeP(S)
    M = MMatrix{3, 3, TC}(ntuple(x -> zero(TC), Val(9)))
    Mof = OffsetArray(M, _basisinds(B), j .+ (-1:1))
    jmm = j-m
    jmmp1 = jmm+1
    jpm = j+m
    jpmp1 = jpm+1
    jp12jp3 = (j+1)*(2j+3)

    Mof[+1, j + 1] = √(T(jmmp1 * (jmmp1 + 1))/2jp12jp3)     * Y[(j+1, m-1)]
    Mof[ 0, j + 1] = -√(T(jmmp1 * jpmp1) / jp12jp3)         * Y[(j+1, m)]
    Mof[-1, j + 1] = √(T(jpmp1 * (jpmp1 + 1)) /2jp12jp3)    * Y[(j+1, m+1)]

    if j > 0
        jjp1 = j*(j+1)
        if m - 1 >= -j
            Mof[+1, j] = -√(T(jpm * jmmp1)/2jjp1) * Y[(j, m-1)]
        end
        Mof[0, j] = m * √(T(1)/jjp1) * Y[(j, m)]
        if m + 1 <= j
            Mof[-1, j] = √(T(jmm * jpmp1)/2jjp1) * Y[(j, m+1)]
        end

        j2jm1 = j*(2j-1)
        if m - 1 >= -(j - 1)
            Mof[+1, j - 1] = √(T(jpm * (jpm - 1))/2j2jm1) * Y[(j-1, m-1)]
        end
        if -(j - 1) <= m <= (j - 1)
            Mof[0, j - 1] = √(T(jmm * jpm)/j2jm1) * Y[(j-1, m)]
        end
        if m + 1 <= j - 1
            Mof[-1, j - 1] = √(T(jmm * (jmm - 1))/2j2jm1) * Y[(j-1, m+1)]
        end
    end

    SMatrix(M)
end
function _vshbasis(::Irreducible, B::Basis, j, m, θ, ϕ, S::VSHCache)
    Y = _vshbasis(Irreducible(), SphericalCovariant(), j, m, θ, ϕ, S)
    C = basisconversionmatrix(SphericalCovariant(), B, θ, ϕ)
    C * Y
end

"""
    vshbasis(Y::AbstractVSH, B::Basis, modes::Union{SphericalHarmonicModes.LM, SphericalHarmonicModes.ML}, θ, ϕ, [S = maximum(SphericalHarmonicModes.l_range(modes))])

Evaluate a set of vector spherical harmonics ``\\mathbf{Y}_{j m}^\\alpha(θ, ϕ)`` for valid values of ``\\alpha``
for all `(j,m)` in `modes`, and return their components in the basis `B`.
A pre-allocated array of scalar spherical harmonics `S` may be passed as the final argument.
"""
function vshbasis(Y::AbstractVSH, B::Basis, modes::Union{ML,LM}, θ, ϕ, S::VSHCache = cache(θ, ϕ, maximum(l_range(modes))))
    v = [vshbasis(Y, B, j, m, θ, ϕ, S) for (j,m) in modes]
    SHArray(v, modes)
end


"""
    genspharm(j, m, θ, ϕ, [S = VectorSphericalHarmonics.cache(θ, ϕ, j)])

Evaluate the diagonal elements of the Phinney-Burridge ([`PB`](@ref)) vector spherical harmonics
in the [`HelicityCovariant`](@ref) basis for one mode `(j,m)`. Returns an array `Y` whose elements
`Y[n]` are given by

```math
Y_{jm}^N(\\theta,\\phi) = \\left[\\mathbf{P}_{j m}^N(θ, ϕ)\\right]^N
```

These components are referred to as "generalized spherical harmonics" by Dahlen & Tromp (1998),
which is the nomenclature used here.

A pre-allocated array of scalar spherical harmonics `S` may be passed as the final argument.
"""
function genspharm(j, m, θ, ϕ, S::VSHCache = cache(θ, ϕ, j))
    Y = vshbasis(PB(), HelicityCovariant(), j, m, θ, ϕ, S)
    M = diag(parent(Y))
    OffsetArray(M, _basisinds(HelicityCovariant()))
end

"""
    genspharm(modes::Union{SphericalHarmonicModes.LM, SphericalHarmonicModes.ML}, θ, ϕ, [S = maximum(SphericalHarmonicModes.l_range(modes))])

Evaluate the diagonal elements of the Phinney-Burridge ([`PB`](@ref)) vector spherical harmonics
in the [`HelicityCovariant`](@ref) basis for all modes `(j,m)` in `modes`.

A pre-allocated array of scalar spherical harmonics `S` may be passed as the final argument.
"""
function genspharm(modes::Union{ML,LM}, θ, ϕ, S::VSHCache = cache(θ, ϕ, maximum(l_range(modes))))
    v = [genspharm(j, m, θ, ϕ, S) for (j,m) in modes]
    SHArray(v, modes)
end

end

"""
    Basis

Abstract supertype of various basis sets that vector spherical harmonics may be decomposed in.
"""
abstract type Basis end
@doc raw"""
    SphericalCovariant <: Basis

The spherical covariant basis ``\chi_\mu`` for ``\mu\in\{-1,0,1\}``
"""
struct SphericalCovariant <: Basis end
@doc raw"""
    Polar <: Basis

The spherical polar basis ``\hat{r}``, ``\hat{\theta}``, ``\hat{\phi}``.
"""
struct Polar <: Basis end
@doc raw"""
    HelicityCovariant  <: Basis

The helicity basis ``\mathbf{e}_\mu`` for ``\mu\in\{-1,0,1\}``
"""
struct HelicityCovariant <: Basis end
@doc raw"""
    Cartesian  <: Basis

The Cartesian basis ``\hat{x}``, ``\hat{y}``, ``\hat{z}``
"""
struct Cartesian <: Basis end

_basisinds(::Union{HelicityCovariant, SphericalCovariant}) = -1:1
_basisinds(::Union{Cartesian, Polar}) = 1:3

#= Define matrices that convert between components C1 and C2 in bases B1 and B2 as
    M * C1 = C2

The corresponding functions are named B1_B2_conversion.
These matrices are conjugates of those that convert between the bases, so conj.(M) * B1 = B2
=#
const HelicitySphericalPolarConversionMatrix = SMatrix{3,3}([
                                        0       1   0
                                        1/√2    0   -1/√2
                                        -im/√2   0    -im/√2
                                        ])

helicity_polar_conversion(θ, ϕ) = HelicitySphericalPolarConversionMatrix
polar_helicity_conversion(θ, ϕ) = helicity_polar_conversion(θ, ϕ)'

function helicity_spherical_conversion(θ, ϕ)
    cisϕ = cis(ϕ)
    sinθ, cosθ = sincos(θ)
    invsqrt2 = 1/√2
    sin²θby2 = (1 - cosθ)/2
    cos²θby2 = (1 + cosθ)/2
    SMatrix{3,3}((
        cos²θby2 * cisϕ,
        -sinθ * invsqrt2,
        sin²θby2 * conj(cisϕ),

        sinθ * invsqrt2 * cisϕ,
        cosθ,
        -sinθ * invsqrt2 * conj(cisϕ),

        sin²θby2 * cisϕ,
        sinθ * invsqrt2,
        cos²θby2 * conj(cisϕ),
        ))
end
spherical_helicity_conversion(θ, ϕ) = helicity_spherical_conversion(θ, ϕ)'

function helicity_cartesian_conversion(θ, ϕ)
    sinθ, cosθ = sincos(θ)
    sinϕ, cosϕ = sincos(ϕ)
    invsqrt2 = 1/√2

    SMatrix{3,3}((
        (cosθ * cosϕ - im * sinϕ) * invsqrt2,
        (cosθ * sinϕ + im * cosϕ) * invsqrt2,
        -sinθ * invsqrt2,

        sinθ * cosϕ,
        sinθ * sinϕ,
        cosθ,

        -(cosθ * cosϕ + im * sinϕ) * invsqrt2,
        -(cosθ * sinϕ - im * cosϕ) * invsqrt2,
        sinθ * invsqrt2,
        ))
end
cartesian_helicity_conversion(θ, ϕ) = helicity_cartesian_conversion(θ, ϕ)'

function spherical_polar_conversion(θ, ϕ)
    invsqrt2 = 1/√2
    invsqrt2cisϕ = invsqrt2 * cis(ϕ)
    sinθ, cosθ = sincos(θ)
    normsinθcisϕ = sinθ * invsqrt2cisϕ
    normcosθcisϕ = cosθ * invsqrt2cisϕ
    SMatrix{3,3}((
        conj(normsinθcisϕ),
        conj(normcosθcisϕ),
        -im * conj(invsqrt2cisϕ),

        cosθ,
        -sinθ,
        0,

        -normsinθcisϕ,
        -normcosθcisϕ,
        -im * invsqrt2cisϕ,
        ))
end
polar_spherical_conversion(θ, ϕ) = spherical_polar_conversion(θ, ϕ)'

const SphericalCartesianConversionMatrix = SMatrix{3,3}((
        1/√2,
        -im/√2,
        0,

        0,
        0,
        1,

        -1/√2,
        -im/√2,
        0,
        ))

spherical_cartesian_conversion(θ, ϕ) = SphericalCartesianConversionMatrix
cartesian_spherical_conversion(θ, ϕ) = spherical_cartesian_conversion(θ, ϕ)'

function cartesian_polar_conversion(θ, ϕ)
    M = basisconversionmatrix(SphericalCovariant(), Polar(), θ, ϕ)
    N = basisconversionmatrix(SphericalCovariant(), Cartesian(), θ, ϕ)
    M * N'
end
polar_cartesian_conversion(θ, ϕ) = cartesian_polar_conversion(θ, ϕ)'

basisconversionmatrix(::T, ::T, θ, ϕ) where {T<:Basis} = I
basisconversionmatrix(::HelicityCovariant, ::SphericalCovariant, θ, ϕ) = helicity_spherical_conversion(θ, ϕ)
basisconversionmatrix(::SphericalCovariant, ::HelicityCovariant, θ, ϕ) = spherical_helicity_conversion(θ, ϕ)
basisconversionmatrix(::HelicityCovariant, ::Polar, θ, ϕ) = helicity_polar_conversion(θ, ϕ)
basisconversionmatrix(::Polar, ::HelicityCovariant, θ, ϕ) = polar_helicity_conversion(θ, ϕ)
basisconversionmatrix(::HelicityCovariant, ::Cartesian, θ, ϕ) = helicity_cartesian_conversion(θ, ϕ)
basisconversionmatrix(::Cartesian, ::HelicityCovariant, θ, ϕ) = cartesian_helicity_conversion(θ, ϕ)

basisconversionmatrix(::SphericalCovariant, ::Polar, θ, ϕ) = spherical_polar_conversion(θ, ϕ)
basisconversionmatrix(::Polar, ::SphericalCovariant, θ, ϕ) = polar_spherical_conversion(θ, ϕ)
basisconversionmatrix(::SphericalCovariant, ::Cartesian, θ, ϕ) = spherical_cartesian_conversion(θ, ϕ)
basisconversionmatrix(::Cartesian, ::SphericalCovariant, θ, ϕ) = cartesian_spherical_conversion(θ, ϕ)

basisconversionmatrix(::Cartesian, ::Polar, θ, ϕ) = cartesian_polar_conversion(θ, ϕ)
basisconversionmatrix(::Polar, ::Cartesian, θ, ϕ) = polar_cartesian_conversion(θ, ϕ)

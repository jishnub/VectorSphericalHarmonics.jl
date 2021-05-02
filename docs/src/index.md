```@meta
CurrentModule = VectorSphericalHarmonics
DocTestSetup = :(using VectorSphericalHarmonics)
```

# VectorSphericalHarmonics.jl

This package lets one compute three variants of vector spherical harmonic (VSH) in different bases following the notation of Varshalovich et al. (1988).

## Vector harmonics

The first (and fundamental) harmonic is the eigenfucntion of the irreducible representatation of total angular momentum ``\mathbf{J} = \mathbf{L}\oplus\mathbf{S}``, with the spin eigenfunctions ``\chi_\mu`` being vectors (corresponding to ``s=1`` and ``\mu\in\{-1,0,1\}``, and satisfying ``S^2\chi_\mu = 2\chi_\mu`` and ``S_z\chi_\mu = \mu\chi_\mu``). The basis spanned by these spin eigenfunctions are referred to as the spherical basis (not to be confused with the spherical polar basis). These vector spherical harmonics may be defined in terms of scalar spherical harmonics ``Y_{L m}\left(\hat{n}\right)`` as

```math
\mathbf{Y}_{J M}^L\left(\hat{n}\right) = \sum_{m \mu} C^{J M}_{L m 1 \mu} Y_{L m}\left(\hat{n}\right) \chi_\mu,
```

where ``C^{J M}_{L m 1 \mu}`` are Clebsch-Gordan coefficients corresponding to the sum of momenta. These vector spherical harmonics satisfy

```math
\begin{aligned}
J^{2}\mathbf{Y}_{JM}^{L}\left(\hat{n}\right) & =J\left(J+1\right)\mathbf{Y}_{JM}^{L}\left(\hat{n}\right),\\
J_{z}\mathbf{Y}_{JM}^{L}\left(\hat{n}\right) & =M\mathbf{Y}_{JM}^{L}\left(\hat{n}\right),\\
L^{2}\mathbf{Y}_{JM}^{L}\left(\hat{n}\right) & =L\left(L+1\right)\mathbf{Y}_{JM}^{L}\left(\hat{n}\right),\\
S^{2}\mathbf{Y}_{JM}^{L}\left(\hat{n}\right) & =2\mathbf{Y}_{JM}^{L}\left(\hat{n}\right).
\end{aligned}
```

In this package we refer to these harmonics as `Irreducible`.

The two other sets of harmonics are linear combinations of the `Irreducible` ones. The first set, referred to as `Hansen` harmonics, are given by

```math
\begin{aligned}
\mathbf{H}_{JM}^{\left(-1\right)}\left(\hat{n}\right) & =-\sqrt{\frac{J+1}{2J+1}}\mathbf{Y}_{JM}^{J+1}\left(\hat{n}\right)+\sqrt{\frac{J}{2J+1}}\mathbf{Y}_{JM}^{J-1}\left(\hat{n}\right),\\
\mathbf{H}_{JM}^{\left(0\right)}\left(\hat{n}\right) & =\mathbf{Y}_{JM}^{J}\left(\hat{n}\right),\\
\mathbf{H}_{JM}^{\left(1\right)}\left(\hat{n}\right) & =\sqrt{\frac{J}{2J+1}}\mathbf{Y}_{JM}^{J+1}\left(\hat{n}\right)+\sqrt{\frac{J+1}{2J+1}}\mathbf{Y}_{JM}^{J-1}\left(\hat{n}\right),
\end{aligned}
```

The Hansen VSH basis is related to scalar spherical harmonics through

```math
\begin{aligned}
\mathbf{H}_{JM}^{\left(-1\right)}\left(\hat{n}\right) & =\mathbf{n}Y_{JM}\left(\hat{n}\right),\\
\mathbf{H}_{JM}^{\left(0\right)}\left(\hat{n}\right) & =\frac{-i}{\sqrt{J\left(J+1\right)}}\left(\mathbf{n}\times\bm{\nabla}_{\Omega}\right)Y_{JM}\left(\hat{n}\right)=\frac{\hat{\mathbf{L}}}{\sqrt{J\left(J+1\right)}}Y_{JM}\left(\hat{n}\right),\\
\mathbf{H}_{JM}^{\left(1\right)}\left(\hat{n}\right) & =\frac{1}{\sqrt{J\left(J+1\right)}}\bm{\nabla}_{\Omega}Y_{JM}\left(\hat{n}\right).
\end{aligned}
```

The third set, referred to as `PB` harmonics following their use by Phinney and Burridge (1973), are related to the `Hansen` harmonics through

```math
\begin{aligned}
\mathbf{P}_{JM}^{+1}\left(\hat{n}\right) & =\frac{1}{\sqrt{2}}\left(\mathbf{H}_{JM}^{\left(1\right)}\left(\hat{n}\right)-\mathbf{H}_{JM}^{\left(0\right)}\left(\hat{n}\right)\right),\\
\mathbf{P}_{JM}^{0}\left(\hat{n}\right) & =\mathbf{H}_{JM}^{\left(-1\right)}\left(\hat{n}\right),\\
\mathbf{P}_{JM}^{-1}\left(\hat{n}\right) & =\frac{1}{\sqrt{2}}\left(\mathbf{H}_{JM}^{\left(1\right)}\left(\hat{n}\right)+\mathbf{H}_{JM}^{\left(0\right)}\left(\hat{n}\right)\right).
\end{aligned}
```

### Properties of VSH

Each set of vector spherical harmonics form a complete, orthonormal basis to decompose 3D vector fields on a sphere.

```math
\int_{0}^{\pi}\int_{0}^{2\pi}\mathbf{Y}_{J^{\prime}M^{\prime}}^{L^{\prime}\dagger}\left(\hat{n}\right)\mathbf{Y}_{JM}^{L}\left(\hat{n}\right)\sin\theta d\theta d\phi=\delta_{J^{\prime}J}\delta_{M^{\prime}M}\delta_{L^{\prime}L}.
```

The various VSH also satisfy some variant of local orthogonality relations. The PB VSH satisfy

```math
\mathbf{P}_{JM}^{\mu\dagger}\left(\hat{n}\right)\cdot\mathbf{P}_{JM}^{\nu}\left(\hat{n}\right)=0,\quad\mu\neq\nu.
```

The Hansen VSH satisfy

```math
\mathbf{H}_{JM}^{\left(\mu\right)}\left(\hat{n}\right)\cdot\mathbf{H}_{JM}^{\left(\nu\right)}\left(\hat{n}\right)=0,\quad\mu\neq\nu.
```

The Irreducible VSH satisfy

```math
\sum_{M}\mathbf{Y}_{JM}^{L^{\prime}\dagger}\left(\hat{n}\right)\cdot\mathbf{Y}_{JM}^{L}\left(\hat{n}\right)=0,\quad L\neq L^{\prime}.
```

In general, the Irreducible VSH satisfies the addition theorem

```math
\sum_{M}\mathbf{Y}_{JM}^{L^{\prime}\dagger}\left(\hat{n}_{1}\right)\cdot\mathbf{Y}_{JM}^{L}\left(\hat{n}_{2}\right)=\delta_{LL^{\prime}}\frac{\left(2J+1\right)}{4\pi}P_{L}\left(\hat{n}_{1}\cdot\hat{n}_{2}\right),
```

where the ``P_L`` are Legendre polynomials. The orthogonality may be seen as a special case of the addition theorem.

## Representation in a basis

There are four different orthonormal bases provided by this package in which the vector spherical harmonics may be represented:

1. [`Cartesian`](@ref) basis (``\hat{x}``, ``\hat{y}`` and ``\hat{z}``)
2. [`Polar`](@ref) basis (``\hat{r}``, ``\hat{\theta}`` and ``\hat{\phi}``)
3. [`SphericalCovariant`](@ref) basis (``\chi_\mu``, that are eigenfunctions of ``S^2`` and ``S_z``, where ``\mu\in\{-1,0,1\}``)
4. [`HelicityCovariant`](@ref) basis (``\mathbf{e}_\mu``, that are eigenfunctions of ``S^2`` and ``\mathbf{S}\cdot\hat{r}``, where ``\mu\in\{-1,0,1\}``)

The vector spherical harmonics ``\mathbf{Y}_{J M}^\alpha`` may be expanded in a basis ``\mathbf{v}_\beta`` as
```math
\mathbf{Y}_{J M}^\alpha = Y_{J M}^{\alpha \beta} \mathbf{v}_\beta,
```
where the components ``Y_{J M}^{\alpha \beta}`` may be expressed as a matrix in the variables ``\alpha`` and ``\beta``. This package evaluates these matrices of coefficients given a harmonic type and a basis set.

The function to evaluate such a matrix is [`vshbasis`](@ref). As an example, to evaluate the components of ``\mathbf{Y}_{10}^1(\pi/3, \pi/3)`` in the spherical polar basis, we may use

```jldoctest
julia> vshbasis(Irreducible(), Polar(), 1, 0, 1, π/3, π/3)
3-element StaticArrays.SVector{3, ComplexF64} with indices SOneTo(3):
  1.6776734621228342e-19 - 1.3874241933809394e-17im
  -4.803764873421188e-18 - 1.0513783175449769e-18im
 -2.4687675144565694e-18 + 0.29920671030107454im
```

In general it's possible to general the entire matrix for one `(J,M)` in one-go, eg.
```jldoctest VSHmatrix
julia> Y = vshbasis(Irreducible(), Polar(), 1, 0, π/3, π/3)
3×3 OffsetArray(::StaticArrays.SMatrix{3, 3, ComplexF64, 9}, 1:3, 0:2) with eltype ComplexF64 with indices 1:3×0:2:
  0.141047+0.0im   1.67767e-19-1.38742e-17im     -0.199471+2.37624e-18im
 -0.244301+0.0im  -4.80376e-18-1.05138e-18im     -0.172747+9.59881e-19im
       0.0+0.0im  -2.46877e-18+0.299207im     -1.16773e-19+1.28643e-17im
```

The matrix elements `Y[α, n]` in this case represent the component ``[\mathbf{Y}_{10}^n(\pi/3, \pi/3)]^\alpha``. For example, the component ``\mathbf{Y}_{10}^1(\pi/3, \pi/3)\cdot\hat{r}`` is given by `Y[1,1]`, the component ``\mathbf{Y}_{10}^1(\pi/3, \pi/3)\cdot\hat{\theta}`` is given by `Y[2,1]` while the the component ``\mathbf{Y}_{10}^0(\pi/3, \pi/3)\cdot\hat{\theta}`` is given by `Y[2,0]`. Note that the VSH superscript always corresponds to the second axis of the matrix, whereas the basis coresponds to the first axis. The indices of the VSH directly correspond those of the matrix, whereas the first axis of the matrix has indices `-1:1` for `SphericalCovariant` and `HelicityCovariant` bases, and `1:3` for `Polar` and `Cartesian` bases.

For the complex bases `SphericalCovariant` and `HelicityCovariant`, the matrix elements represent the contravariant components of the vector harmonics in the respective basis.

The PB basis is of particular importance, as it is locally diagonal in the helicity basis.

```jldoctest
julia> vshbasis(PB(), HelicityCovariant(), 1, 0, π/3, π/3)
3×3 OffsetArray(::LinearAlgebra.Diagonal{ComplexF64, StaticArrays.SVector{3, ComplexF64}}, -1:1, -1:1) with eltype ComplexF64 with indices -1:1×-1:1:
 -0.299207-1.79592e-18im           ⋅               ⋅
           ⋅              0.244301+0.0im           ⋅
           ⋅                       ⋅      0.299207-4.35925e-18im
```

The vectors harmonics themselves, therefore, are orthogonal at each point. One may obtain the diagonal elements through the function [`genspharm`](@ref).

```jldoctest
julia> genspharm(1, 0, π/3, π/3)
3-element OffsetArray(::StaticArrays.SVector{3, ComplexF64}, -1:1) with eltype ComplexF64 with indices -1:1:
 -0.2992067103010745 - 1.7959178942769708e-18im
 0.24430125595146002 + 0.0im
  0.2992067103010745 - 4.359250168826542e-18im
```

The diagonal elements of the PB VSH basis are related to the Wigner d-matrix through

```math
\left[\mathbf{P}_{JM}^{\alpha}\left(\theta,\phi\right)\right]^{\alpha}=\sqrt{\frac{2J+1}{4\pi}}d_{M\alpha}^{J}\left(\theta\right)\exp\left(iM\phi\right),
```

and the elements for ``\alpha=0`` are scalar spherical harmonics

```math
\left[\mathbf{Y}_{JM}^{0}\left(\theta,\phi\right)\right]^{0}=Y_{JM}\left(\theta,\phi\right).
```

One may also compute the matrices for a range of modes using the iterators provided by [`SphericalHarmonicModes.jl`](https://github.com/jishnub/SphericalHarmonicModes.jl). As an example, we may evaluate the matrices for all `M` for `J = 1:3` as
```jldoctest VSHmodes
julia> using SphericalHarmonicModes

julia> Y = vshbasis(Irreducible(), Polar(), ML(1:3), π/3, π/3);
```

This returns a vector that may be indexed using `(j,m)` to obtain the individual component matrices as

```jldoctest VSHmodes
julia> Y[(1,0)]
3×3 OffsetArray(::StaticArrays.SMatrix{3, 3, ComplexF64, 9}, 1:3, 0:2) with eltype ComplexF64 with indices 1:3×0:2:
  0.141047+0.0im   1.67767e-19-1.38742e-17im     -0.199471+2.37624e-18im
 -0.244301+0.0im  -4.80376e-18-1.05138e-18im     -0.172747+9.59881e-19im
       0.0+0.0im  -2.46877e-18+0.299207im     -1.16773e-19+1.28643e-17im
```

## Pre-allocation

For performance reasons, it is prudent to pre-allocate a set of scalar spherical harmonics that are used to compute the vector harmonics. This may be carried out using the function [`cache`](@ref) as

```jldoctest cache
julia> S = VectorSphericalHarmonics.cache(π/3, π/3, 30);
```

This cache may be passed to `vshbasis` as the last argument as

```jldoctest cache
julia> vshbasis(PB(), HelicityCovariant(), 1, 0, π/3, π/3, S)
3×3 OffsetArray(::LinearAlgebra.Diagonal{ComplexF64, StaticArrays.SVector{3, ComplexF64}}, -1:1, -1:1) with eltype ComplexF64 with indices -1:1×-1:1:
 -0.299207-1.79592e-18im           ⋅               ⋅
           ⋅              0.244301+0.0im           ⋅
           ⋅                       ⋅      0.299207-4.35925e-18im
```

This will help performance as repeated evaluations may be avoided.

```julia
julia> @btime vshbasis(PB(), HelicityCovariant(), 1, 0, π/3, π/3, $S);
  201.503 ns (0 allocations: 0 bytes)
```

# Poles

There are special methods defined to compute vector spherical harmonics at the poles. One may use the `NorthPole` and `SouthPole` types defined in [`SphericalHarmonics.jl`](https://github.com/jishnub/SphericalHarmonics.jl/) to evaluate these.

```jldoctest
julia> using SphericalHarmonics: NorthPole

julia> vshbasis(Irreducible(), Polar(), 1, 0, NorthPole(), 0)
3×3 OffsetArray(::StaticArrays.SMatrix{3, 3, ComplexF64, 9}, 1:3, 0:2) with eltype ComplexF64 with indices 1:3×0:2:
 0.282095+0.0im  0.0+0.0im  -0.398942+0.0im
     -0.0+0.0im  0.0+0.0im        0.0+0.0im
      0.0+0.0im  0.0+0.0im        0.0+0.0im
```

The vector spherical harmonics are non-zero at the poles only for ``m=0,\pm 1``.

# Rotation of coordinates

We assume that a point ``\hat{n}`` has the coordinates ``(\theta_1,\phi_1)`` in the frame ``S_1`` and ``(\theta_2,\phi_2)`` in the frame ``S_2``, where the two frames are related by a rotation ``S_2 = R_z\left(\phi_1\right)R_y\left(\theta_1-\theta_2\right)R_z\left(-\phi_2\right)`` about a set of axes fixed in space. Under this rotation, vector spherical harmonics at the two points are related by

```math
\mathbf{Y}_{JM^{\prime}}^{\alpha}\left(\theta_{2},\phi_{2}\right)=\sum_{M}D_{MM^{\prime}}^{J}\left(\phi_{1},\theta_{1}-\theta_{2},-\phi_{2}\right)\mathbf{Y}_{JM}^{\alpha}\left(\theta_{1},\phi_{1}\right),
```

where ``D_{MM^{\prime}}^{J}`` are elements of the Wigner D-matrix. This relation holds for all the harmonics defined here. An important point to note here is that the new vector spherical harmonics ``\mathbf{Y}_{JM^{\prime}}^{\alpha}\left(\theta_{2},\phi_{2}\right)`` are computed about the rotated set of axes ``S_2``, whereas ``\mathbf{Y}_{JM}^{\alpha}\left(\theta_{1},\phi_{1}\right)`` are evaluated about ``S_1``. A similar relation may be obtained for the `Polar` and `HelicityCovariant` components of vector spherical harmonics.

!!! note
    This relation does not hold for the `Cartesian` and `SphericalCovariant` components of the vector harmonics, as the axes of ``S_1`` and ``S_2`` are oriented differently in space.

# Index

```@autodocs
Modules = [VectorSphericalHarmonics]
```

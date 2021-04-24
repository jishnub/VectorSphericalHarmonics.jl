var documenterSearchIndex = {"docs":
[{"location":"","page":"Reference","title":"Reference","text":"CurrentModule = VectorSphericalHarmonics\nDocTestSetup = :(using VectorSphericalHarmonics)","category":"page"},{"location":"#VectorSphericalHarmonics.jl","page":"Reference","title":"VectorSphericalHarmonics.jl","text":"","category":"section"},{"location":"","page":"Reference","title":"Reference","text":"This package lets one compute three variants of vector spherical harmonic (VSH) in different bases following the notation of Varshalovich et al. (1988).","category":"page"},{"location":"#Vector-harmonics","page":"Reference","title":"Vector harmonics","text":"","category":"section"},{"location":"","page":"Reference","title":"Reference","text":"The first (and fundamental) harmonic is the eigenfucntion of the irreducible representatation of total angular momentum mathbfJ = mathbfLoplusmathbfS, with the spin eigenfunctions chi_mu being vectors (corresponding to s=1 and muin-101, and satisfying S^2chi_mu = 2chi_mu and S_zchi_mu = muchi_mu). The basis spanned by these spin eigenfunctions are referred to as the spherical basis (not to be confused with the spherical polar basis). These vector spherical harmonics may be defined in terms of scalar spherical harmonics Y_L mleft(hatnright) as","category":"page"},{"location":"","page":"Reference","title":"Reference","text":"mathbfY_J M^Lleft(hatnright) = sum_m mu C^J M_L m 1 mu Y_L mleft(hatnright) chi_mu","category":"page"},{"location":"","page":"Reference","title":"Reference","text":"where C^J M_L m 1 mu are Clebsch-Gordan coefficients corresponding to the sum of momenta. These vector spherical harmonics satisfy","category":"page"},{"location":"","page":"Reference","title":"Reference","text":"beginaligned\nJ^2mathbfY_JM^Lleft(hatnright)  =Jleft(J+1right)mathbfY_JM^Lleft(hatnright)\nJ_zmathbfY_JM^Lleft(hatnright)  =MmathbfY_JM^Lleft(hatnright)\nL^2mathbfY_JM^Lleft(hatnright)  =Lleft(L+1right)mathbfY_JM^Lleft(hatnright)\nS^2mathbfY_JM^Lleft(hatnright)  =2mathbfY_JM^Lleft(hatnright)\nendaligned","category":"page"},{"location":"","page":"Reference","title":"Reference","text":"In this package we refer to these harmonics as Irreducible.","category":"page"},{"location":"","page":"Reference","title":"Reference","text":"The two other sets of harmonics are linear combinations of the Irreducible ones. The first set, referred to as Hansen harmonics, are given by","category":"page"},{"location":"","page":"Reference","title":"Reference","text":"beginaligned\nmathbfH_JM^left(-1right)left(hatnright)  =-sqrtfracJ+12J+1mathbfY_JM^J+1left(hatnright)+sqrtfracJ2J+1mathbfY_JM^J-1left(hatnright)\nmathbfH_JM^left(0right)left(hatnright)  =mathbfY_JM^Jleft(hatnright)\nmathbfH_JM^left(1right)left(hatnright)  =sqrtfracJ2J+1mathbfY_JM^J+1left(hatnright)+sqrtfracJ+12J+1mathbfY_JM^J-1left(hatnright)\nendaligned","category":"page"},{"location":"","page":"Reference","title":"Reference","text":"The third set, referred to as PB harmonics following their use by Phinney and Burridge (1973), are related to the Hansen harmonics through","category":"page"},{"location":"","page":"Reference","title":"Reference","text":"beginaligned\nmathbfP_JM^+1left(hatnright)  =frac1sqrt2left(mathbfH_JM^left(1right)left(hatnright)-mathbfH_JM^left(0right)left(hatnright)right)\nmathbfP_JM^0left(hatnright)  =mathbfH_JM^left(-1right)left(hatnright)\nmathbfP_JM^-1left(hatnright)  =frac1sqrt2left(mathbfH_JM^left(1right)left(hatnright)+mathbfH_JM^left(0right)left(hatnright)right)\nendaligned","category":"page"},{"location":"","page":"Reference","title":"Reference","text":"Each set of vector spherical harmonics form a complete, orthonormal basis to decompose 3D vector fields on a sphere.","category":"page"},{"location":"#Representation-in-a-basis","page":"Reference","title":"Representation in a basis","text":"","category":"section"},{"location":"","page":"Reference","title":"Reference","text":"There are four different bases of vectors in which the vector spherical harmonics may be represented:","category":"page"},{"location":"","page":"Reference","title":"Reference","text":"Cartesian basis (hatx, haty and hatz)\nSpherical polar basis (hatr, hattheta and hatphi)\nSpherical (chi_mu, that are eigenfunctions of S^2 and S_z, where muin-101)\nHelicity (mathbfe_mu, that are eigenfunctions of S^2 and mathbfScdothatr, where muin-101)","category":"page"},{"location":"","page":"Reference","title":"Reference","text":"The vector spherical harmonics mathbfY_J M^alpha may be expanded in a basis mathbfv_beta as","category":"page"},{"location":"","page":"Reference","title":"Reference","text":"mathbfY_J M^alpha = Y_J M^alpha beta mathbfv_beta","category":"page"},{"location":"","page":"Reference","title":"Reference","text":"where the components Y_J M^alpha beta may be expressed as a matrix in the variables alpha and beta. This package evaluates these matrices of coefficients given a harmonic type and a basis set.","category":"page"},{"location":"","page":"Reference","title":"Reference","text":"The function to evaluate such a matrix is vshbasis. As an example, to evaluate the components of Y_10^1(pi3 pi3) in the spherical polar basis, we may use","category":"page"},{"location":"","page":"Reference","title":"Reference","text":"julia> vshbasis(Irreducible(), Polar(), 1, 0, 1, π/3, π/3)\n3-element StaticArrays.SVector{3, ComplexF64} with indices SOneTo(3):\n    -0.19947114020071638 + 2.3762446998036102e-18im\n    -0.17274707473566778 + 9.59881327122146e-19im\n -1.1677330786160424e-19 + 1.2864343079066462e-17im","category":"page"},{"location":"","page":"Reference","title":"Reference","text":"In general it's possible to general the entire matrix for one (J,M) in one-go, eg.","category":"page"},{"location":"","page":"Reference","title":"Reference","text":"julia> Y = vshbasis(Irreducible(), Polar(), 1, 0, π/3, π/3)\n3×3 OffsetArray(::StaticArrays.SMatrix{3, 3, ComplexF64, 9}, 1:3, 0:2) with eltype ComplexF64 with indices 1:3×0:2:\n  0.141047+0.0im   1.67767e-19-1.38742e-17im     -0.199471+2.37624e-18im\n -0.244301+0.0im  -4.80376e-18-1.05138e-18im     -0.172747+9.59881e-19im\n       0.0+0.0im  -2.46877e-18+0.299207im     -1.16773e-19+1.28643e-17im","category":"page"},{"location":"","page":"Reference","title":"Reference","text":"Each column of this matrix represents one vector, and each row represents the projection of the VSH along the spherical polar unit vectors (hatr, hattheta, hatphi in order). For example, the component Y_10^1(pi3 pi3)cdothatr is given by Y[1,1], the component Y_10^1(pi3 pi3)cdothattheta is given by Y[2,1] while the the component Y_10^0(pi3 pi3)cdothattheta is given by Y[2,0]. Note that the VSH superscript always corresponds to the second axis of the matrix, whereas the basis coresponds to the first axis. The indices of the VSH directly correspond those of the matrix, whereas the first axis of the matrix has indices 0:2 for Spherical and Helicity bases, and 1:3 for Spherical and Cartesian bases.","category":"page"},{"location":"","page":"Reference","title":"Reference","text":"The PB basis is of particular importance, as it is locally diagonal in the helicity basis.","category":"page"},{"location":"","page":"Reference","title":"Reference","text":"julia> vshbasis(PB(), HelicityCovariant(), 1, 0, π/3, π/3)\n3×3 OffsetArray(::LinearAlgebra.Diagonal{ComplexF64, StaticArrays.SVector{3, ComplexF64}}, -1:1, -1:1) with eltype ComplexF64 with indices -1:1×-1:1:\n -0.299207-1.79592e-18im           ⋅               ⋅\n           ⋅              0.244301+0.0im           ⋅\n           ⋅                       ⋅      0.299207-4.35925e-18im","category":"page"},{"location":"","page":"Reference","title":"Reference","text":"The vectors themselves, therefore, are orthogonal at each point. Such a relation does not hold for the other harmonics.","category":"page"},{"location":"","page":"Reference","title":"Reference","text":"One may also compute the matrices for a range of modes using the iterators provided by SphericalHarmonicModes.jl. As an example, we may evaluate the matrices for all M for J = 1:3 as","category":"page"},{"location":"","page":"Reference","title":"Reference","text":"julia> using SphericalHarmonicModes\n\njulia> Y = vshbasis(Irreducible(), Polar(), ML(1:3), π/3, π/3);","category":"page"},{"location":"","page":"Reference","title":"Reference","text":"This returns a vector that may be indexed using (j,m) to obtain the individual component matrices as","category":"page"},{"location":"","page":"Reference","title":"Reference","text":"julia> Y[(1,0)]\n3×3 OffsetArray(::StaticArrays.SMatrix{3, 3, ComplexF64, 9}, 1:3, 0:2) with eltype ComplexF64 with indices 1:3×0:2:\n  0.141047+0.0im   1.67767e-19-1.38742e-17im     -0.199471+2.37624e-18im\n -0.244301+0.0im  -4.80376e-18-1.05138e-18im     -0.172747+9.59881e-19im\n       0.0+0.0im  -2.46877e-18+0.299207im     -1.16773e-19+1.28643e-17im","category":"page"},{"location":"#Pre-allocation","page":"Reference","title":"Pre-allocation","text":"","category":"section"},{"location":"","page":"Reference","title":"Reference","text":"For performance reasons, it is prudent to pre-allocate a set of scalar spherical harmonics that are used to compute the vector harmonics. This may be carried out using the function cache as","category":"page"},{"location":"","page":"Reference","title":"Reference","text":"julia> S = VectorSphericalHarmonics.cache(π/3, π/3, 30);","category":"page"},{"location":"","page":"Reference","title":"Reference","text":"This cache may be passed to vshbasis as the last argument as","category":"page"},{"location":"","page":"Reference","title":"Reference","text":"julia> vshbasis(PB(), HelicityCovariant(), 1, 0, π/3, π/3, S)\n3×3 OffsetArray(::LinearAlgebra.Diagonal{ComplexF64, StaticArrays.SVector{3, ComplexF64}}, -1:1, -1:1) with eltype ComplexF64 with indices -1:1×-1:1:\n -0.299207-1.79592e-18im           ⋅               ⋅\n           ⋅              0.244301+0.0im           ⋅\n           ⋅                       ⋅      0.299207-4.35925e-18im","category":"page"},{"location":"","page":"Reference","title":"Reference","text":"This will help performance as repeated evaluations may be avoided.","category":"page"},{"location":"","page":"Reference","title":"Reference","text":"julia> @btime vshbasis(PB(), HelicityCovariant(), 1, 0, π/3, π/3, $S);\n  201.503 ns (0 allocations: 0 bytes)","category":"page"},{"location":"#Poles","page":"Reference","title":"Poles","text":"","category":"section"},{"location":"","page":"Reference","title":"Reference","text":"There are special methods defined to compute vector spherical harmonics at the poles. One may use the NorthPole and SouthPole types defined in SphericalHarmonics.jl to evaluate these.","category":"page"},{"location":"","page":"Reference","title":"Reference","text":"julia> using SphericalHarmonics: NorthPole\n\njulia> vshbasis(Irreducible(), Polar(), 1, 0, NorthPole(), 0)\n3×3 OffsetArray(::StaticArrays.SMatrix{3, 3, ComplexF64, 9}, 1:3, 0:2) with eltype ComplexF64 with indices 1:3×0:2:\n 0.282095+0.0im  0.0+0.0im  -0.398942+0.0im\n     -0.0+0.0im  0.0+0.0im        0.0+0.0im\n      0.0+0.0im  0.0+0.0im        0.0+0.0im","category":"page"},{"location":"","page":"Reference","title":"Reference","text":"The vector spherical harmonics are non-zero at the poles only for m=0pm 1.","category":"page"},{"location":"#Index","page":"Reference","title":"Index","text":"","category":"section"},{"location":"","page":"Reference","title":"Reference","text":"Modules = [VectorSphericalHarmonics]","category":"page"},{"location":"#VectorSphericalHarmonics.AbstractVSH","page":"Reference","title":"VectorSphericalHarmonics.AbstractVSH","text":"AbstractVSH\n\nAbstract supertype of vector spherical harmonics\n\n\n\n\n\n","category":"type"},{"location":"#VectorSphericalHarmonics.Basis","page":"Reference","title":"VectorSphericalHarmonics.Basis","text":"Basis\n\nAbstract supertype of various basis sets that vector spherical harmonics may be decomposed in.\n\n\n\n\n\n","category":"type"},{"location":"#VectorSphericalHarmonics.Cartesian","page":"Reference","title":"VectorSphericalHarmonics.Cartesian","text":"Cartesian  <: Basis\n\nThe Cartesian basis hatx, haty, hatz\n\n\n\n\n\n","category":"type"},{"location":"#VectorSphericalHarmonics.Hansen","page":"Reference","title":"VectorSphericalHarmonics.Hansen","text":"Hansen <: AbstractVSH\n\nHansen vector spherical harmonics mathbfH_J M^(lambda)left(hatnright) that are related to the Irreducible harmonics through\n\nbeginaligned\nmathbfH_JM^left(-1right)left(hatnright)  =-sqrtfracJ+12J+1mathbfY_JM^J+1left(hatnright)+sqrtfracJ2J+1mathbfY_JM^J-1left(hatnright)\nmathbfH_JM^left(0right)left(hatnright)  =mathbfY_JM^Jleft(hatnright)\nmathbfH_JM^left(1right)left(hatnright)  =sqrtfracJ2J+1mathbfY_JM^J+1left(hatnright)+sqrtfracJ+12J+1mathbfY_JM^J-1left(hatnright)\nendaligned\n\n\n\n\n\n","category":"type"},{"location":"#VectorSphericalHarmonics.HelicityCovariant","page":"Reference","title":"VectorSphericalHarmonics.HelicityCovariant","text":"HelicityCovariant  <: Basis\n\nThe helicity basis mathbfe_mu for muin-101\n\n\n\n\n\n","category":"type"},{"location":"#VectorSphericalHarmonics.Irreducible","page":"Reference","title":"VectorSphericalHarmonics.Irreducible","text":"Irreducible <: AbstractVSH\n\nVector spherical harmonics that are eigenfunctions of irreducible representationa of total, orbital as well as spin angular momenta. They may be constructed by coupling scalar spherical harmonics Y_L mleft(hatnright) with the spherical basis vectors chi_mu as\n\nmathbfY_J M^Lleft(hatnright) = sum_m mu C^J M_L m 1 mu Y_L mleft(hatnright) chi_mu\n\n\n\n\n\n","category":"type"},{"location":"#VectorSphericalHarmonics.PB","page":"Reference","title":"VectorSphericalHarmonics.PB","text":"PB <: AbstractVSH\n\nPhinney-Burridge vector spherical harmonics mathbfP_J M^gammaleft(hatnright) that are related to the Hansen harmonics through\n\nbeginaligned\nmathbfP_JM^+1left(hatnright)  =frac1sqrt2left(mathbfH_JM^left(1right)left(hatnright)-mathbfH_JM^left(0right)left(hatnright)right)\nmathbfP_JM^0left(hatnright)  =mathbfH_JM^left(-1right)left(hatnright)\nmathbfP_JM^-1left(hatnright)  =frac1sqrt2left(mathbfH_JM^left(1right)left(hatnright)+mathbfH_JM^left(0right)left(hatnright)right)\nendaligned\n\n\n\n\n\n","category":"type"},{"location":"#VectorSphericalHarmonics.Polar","page":"Reference","title":"VectorSphericalHarmonics.Polar","text":"Polar <: Basis\n\nThe spherical polar basis hatr, hattheta, hatphi.\n\n\n\n\n\n","category":"type"},{"location":"#VectorSphericalHarmonics.SphericalCovariant","page":"Reference","title":"VectorSphericalHarmonics.SphericalCovariant","text":"SphericalCovariant <: Basis\n\nThe spherical covariant basis chi_mu for muin-101\n\n\n\n\n\n","category":"type"},{"location":"#VectorSphericalHarmonics.cache!-Tuple{VectorSphericalHarmonics.VSHCache, Any, Any, Vararg{Any, N} where N}","page":"Reference","title":"VectorSphericalHarmonics.cache!","text":"cache!(S, θ, ϕ, [jmax])\n\nUpdate a pre-allocated set of scalar spherical harmonics Y_j m(thetaphi). The cutoff jmax corresponds to the maximum degree of vector harmonics that we seek to evaluate. If this is not provided, the cutoff angular degree in S will remain unchanged.\n\n\n\n\n\n","category":"method"},{"location":"#VectorSphericalHarmonics.cache-Tuple{Any, Any, Any}","page":"Reference","title":"VectorSphericalHarmonics.cache","text":"cache([T::Type = Float64], θ, ϕ, jmax)\n\nPre-allocate a set of scalar spherical harmonics Y_j m(thetaphi) that go into evaluating the vector harmonics. The cutoff jmax corresponds to the maximum degree of vector harmonics that we seek to evaluate.\n\nThe type T sets the precision used to evaluate the harmonics. The spherical harmonics evaluated will be of type Complex{T}.\n\n\n\n\n\n","category":"method"},{"location":"#VectorSphericalHarmonics.vshbasis","page":"Reference","title":"VectorSphericalHarmonics.vshbasis","text":"vshbasis(Y::AbstractVSH, B::Basis, j::Integer, m::Integer, n::Integer, θ, ϕ, [S = VectorSphericalHarmonics.cache(θ, ϕ, j)])\n\nEvaluate the components of the vector spherical harmonics Y_j m^n(θ ϕ) in the basis B. A pre-allocated array of scalar spherical harmonics S may be passed as the final argument.\n\n\n\n\n\n","category":"function"},{"location":"#VectorSphericalHarmonics.vshbasis-2","page":"Reference","title":"VectorSphericalHarmonics.vshbasis","text":"vshbasis(Y::AbstractVSH, B::Basis, modes::Union{SphericalHarmonicModes.LM, SphericalHarmonicModes.ML}, θ, ϕ, [S = maximum(SphericalHarmonicModes.l_range(modes))])\n\nEvaluate a set of vector spherical harmonics Y_j m^alpha(θ ϕ) for valid values of alpha for all (j,m) in modes, and return their components in the basis B. A pre-allocated array of scalar spherical harmonics S may be passed as the final argument.\n\n\n\n\n\n","category":"function"},{"location":"#VectorSphericalHarmonics.vshbasis-3","page":"Reference","title":"VectorSphericalHarmonics.vshbasis","text":"vshbasis(Y::AbstractVSH, B::Basis, j::Integer, m::Integer, θ, ϕ, [S = VectorSphericalHarmonics.cache(θ, ϕ, j)])\n\nEvaluate a set of vector spherical harmonics Y_j m^alpha(θ ϕ) for valid values of alpha, and return their components in the basis B. A pre-allocated array of scalar spherical harmonics S may be passed as the final argument.\n\n\n\n\n\n","category":"function"}]
}

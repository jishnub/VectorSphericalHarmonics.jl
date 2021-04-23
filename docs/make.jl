using Documenter
using VectorSphericalHarmonics

DocMeta.setdocmeta!(VectorSphericalHarmonics, :DocTestSetup, :(using VectorSphericalHarmonics); recursive=true)

makedocs(;
    modules=[VectorSphericalHarmonics],
    authors="Jishnu Bhattacharya",
    repo="https://github.com/jishnub/VectorSphericalHarmonics.jl/blob/{commit}{path}#L{line}",
    sitename="VectorSphericalHarmonics.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://jishnub.github.io/VectorSphericalHarmonics.jl",
        assets=String[],
    ),
    pages=[
        "Reference" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/jishnub/VectorSphericalHarmonics.jl",
)

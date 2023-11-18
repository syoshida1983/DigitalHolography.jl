using DigitalHolography
using Documenter

DocMeta.setdocmeta!(DigitalHolography, :DocTestSetup, :(using DigitalHolography); recursive=true)

makedocs(;
    modules=[DigitalHolography],
    authors="Shuhei Yoshida <yshuhei@ele.kindai.ac.jp> and contributors",
    repo="https://github.com/syoshida1983/DigitalHolography.jl/blob/{commit}{path}#{line}",
    sitename="DigitalHolography.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://syoshida1983.github.io/DigitalHolography.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/syoshida1983/DigitalHolography.jl",
    devbranch="master",
)

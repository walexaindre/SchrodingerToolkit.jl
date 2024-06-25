using SchrodingerToolkit
using Documenter

DocMeta.setdocmeta!(SchrodingerToolkit, :DocTestSetup, :(using SchrodingerToolkit); recursive=true)

makedocs(;
    modules=[SchrodingerToolkit],
    authors="walexaindre <walexaindre@hotmail.com
> and contributors",
    sitename="SchrodingerToolkit.jl",
    format=Documenter.HTML(;
        canonical="https://walexaindre.github.io/SchrodingerToolkit.jl",
        edit_link="dev",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/walexaindre/SchrodingerToolkit.jl",
    devbranch="dev",
)

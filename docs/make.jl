using ControlCenter
using Documenter

DocMeta.setdocmeta!(ControlCenter, :DocTestSetup, :(using ControlCenter); recursive=true)

makedocs(;
    modules=[ControlCenter],
    authors="Malte Krug <malte.krug@uni-ulm.de> and contributors",
    sitename="ControlCenter.jl",
    format=Documenter.HTML(;
        canonical="https://krugx.github.io/ControlCenter.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/krugx/ControlCenter.jl",
    devbranch="main",
)

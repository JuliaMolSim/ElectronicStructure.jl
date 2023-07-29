using ElectronicStructure
using Documenter

DocMeta.setdocmeta!(ElectronicStructure, :DocTestSetup, :(using ElectronicStructure); recursive=true)

makedocs(;
    modules=[ElectronicStructure],
    authors="Michael F. Herbst <info@michael-herbst.com> and contributors",
    repo="https://github.com/mfherbst/ElectronicStructure.jl/blob/{commit}{path}#{line}",
    sitename="ElectronicStructure.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://mfherbst.github.io/ElectronicStructure.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/mfherbst/ElectronicStructure.jl",
    devbranch="master",
)

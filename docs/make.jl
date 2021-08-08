using Glissa
using Documenter

DocMeta.setdocmeta!(Glissa, :DocTestSetup, :(using Glissa); recursive=true)

makedocs(;
    modules=[Glissa],
    authors="Andrew Dolgert <adolgert@uw.edu>",
    repo="https://github.com/adolgert/Glissa.jl/blob/{commit}{path}#{line}",
    sitename="Glissa.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://adolgert.github.io/Glissa.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/adolgert/Glissa.jl",
    devbranch="main",
)

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
        "Notation" => "notation.md",
        "B-splines" => "bsplines.md",
        "Testing B-splines" => "testspaces.md",
        "Reference" => "reference.md",
        "Project Plan" => "projectplan.md",
        "Test Plan" => "testplan.md"
    ],
)

deploydocs(;
    repo="github.com/adolgert/Glissa.jl",
    devbranch="main",
)

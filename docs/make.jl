using Documenter
using BSE_SC

makedocs(
    modules = [BSE_SC]
    sitename = "Selfconsistent Asymptotics Improvement",
    authors="Julian Stobbe <Atomtomate@gmx.de> and contributors",
    repo="https://github.com/Atomtomate/BSE_SC.jl/blob/{commit}{path}#L{line}",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", nothing) == "true",
        canonical="https://Atomtomate.github.io/BSE_SC.jl",
        assets=String[],
    ),
    pages = ["Home" => "index.md"]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    branch="gh-pages",
    devbranch = "master",
    devurl = "stable",
    repo = "github.com/Atomtomate/BSE_SC.jl.git"
)

using Documenter
using e3nn

mathengine = MathJax3()
prettyurls = get(ENV, "CI", nothing) == "true"

makedocs(;
    sitename = "e3nn.jl",
    doctest = false,
    clean = true,
    format = Documenter.HTML(;
        mathengine,
        prettyurls,
        # assets = assets,
        size_threshold = nothing,
    ),
    pages = ["Home" => "index.md", "API Reference" => ["o3" => "api/o3.md"]],
)

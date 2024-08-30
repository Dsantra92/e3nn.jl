using Documenter, e3nn

mathengine = MathJax3()
prettyurls = get(ENV, "CI", nothing) == "true"

makedocs(
    modules = [e3nn],
    doctest = false,
    clean = true,
    sitename = "e3nn.jl",
    format = Documenter.HTML(;
        mathengine,
        prettyurls,
        assets = ["assets/favicon.ico"],
        collapselevel = 3
    ),
    pages = ["Home" => "index.md", "API Reference" => ["o3" => "api/o3.md"]]
)

deploydocs(repo = "github.com/Dsantra92/e3nn.jl.git")

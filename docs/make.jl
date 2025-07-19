using Documenter, E3NN

mathengine = MathJax3()
prettyurls = get(ENV, "CI", nothing) == "true"

makedocs(
    modules = [E3NN],
    doctest = false,
    clean = true,
    sitename = "E3NN.jl",
    format = Documenter.HTML(;
        mathengine,
        prettyurls,
        assets = ["assets/favicon.ico"],
        collapselevel = 3
    ),
    pages = ["Home" => "index.md", "Irreps" => "irreps.md",
        "Examples" => "examples.md",
        "API Reference" => ["Irreps" => "api/irreps.md"]]
)

deploydocs(repo = "github.com/Dsantra92/e3nn.jl.git")

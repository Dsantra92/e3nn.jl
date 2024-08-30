using Documenter, e3nn

mathengine = MathJax3()

makedocs(
    modules = [e3nn],
    doctest = false,
    clean = true,
    sitename="e3nn.jl",
    format = Documenter.HTML(
        canonical = "https://dsantra92.github.io/e3nn.jl/stable/",
        mathengine,
        assets = ["assets/favicon.ico"],
        prettyurls = get(ENV, "CI", nothing) == "true",
        collapselevel=3,
    ),
        pages = ["Home" => "Index.md",
             "API Reference" => [
                 "o3" => "api/o3.md",
             ],
    ]
    )

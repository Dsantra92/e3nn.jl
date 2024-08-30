using Documenter, e3nn

makedocs(
    modules = [e3nn],
    doctest = false,
    clean = false,
    sitename="e3nn.jl",
    format = Documenter.HTML(
        canonical = "https://dsantra92.github.io/e3nn.jl/stable/",
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

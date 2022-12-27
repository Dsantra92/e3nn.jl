using Documenter, e3nn

makedocs(
    modules = [MLDatasets],
    doctest = true,
    clean = false,
    sitename="e3nn.jl",
    format = Documenter.HTML(
        canonical = "https://dsantra92.github.io/e3nn.jl/stable/",
        assets = ["assets/favicon.ico"],
        prettyurls = get(ENV, "CI", nothing) == "true",
        collapselevel=3,
    ),
    )
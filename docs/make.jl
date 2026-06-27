using Pkg

Pkg.activate(@__DIR__)
Pkg.develop(PackageSpec(path=joinpath(@__DIR__, "..")))
Pkg.instantiate()

using Documenter
using UnsteadyKineticRotorDynamics

DocMeta.setdocmeta!(UnsteadyKineticRotorDynamics, :DocTestSetup, :(using CCBlade, UnsteadyKineticRotorDynamics);
    recursive=true)

is_ci = get(ENV, "CI", "false") == "true"

doc_kwargs = Dict(
    :modules => [UnsteadyKineticRotorDynamics],
    :sitename => "UnsteadyKineticRotorDynamics.jl",
    :format => Documenter.HTML(
        prettyurls = is_ci,
        edit_link = "master",
        repolink = "https://github.com/sandialabs/UnsteadyKineticRotorDynamics.jl",
    ),
    :pages => [
        "Home" => "index.md",
        "Quickstart" => "quickstart.md",
        "Theory" => "theory.md",
        "API" => "api.md",
    ],
)

if is_ci
    doc_kwargs[:repo] = "github.com/sandialabs/UnsteadyKineticRotorDynamics.jl"
else
    doc_kwargs[:remotes] = nothing
end

makedocs(; doc_kwargs...)

if is_ci
    deploydocs(
        repo="github.com/sandialabs/UnsteadyKineticRotorDynamics.jl.git",
        devbranch="master",
    )
end

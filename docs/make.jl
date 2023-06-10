using Documenter, Literate, CairoMakie

using GeophysicalFlows

#####
##### Generate literated examples
#####

const EXAMPLES_DIR = joinpath(@__DIR__, "..", "examples")
const OUTPUT_DIR   = joinpath(@__DIR__, "src/literated")

examples = [
  "twodnavierstokes_decaying.jl",
  "twodnavierstokes_stochasticforcing.jl",
  "twodnavierstokes_stochasticforcing_budgets.jl",
  "singlelayerqg_betadecay.jl",
  "singlelayerqg_betaforced.jl",
  "singlelayerqg_decaying_topography.jl",
  "singlelayerqg_decaying_barotropic_equivalentbarotropic.jl",
  "barotropicqgql_betaforced.jl",
  "multilayerqg_2layer.jl",
  "surfaceqg_decaying.jl",
]


for example in examples
  withenv("GITHUB_REPOSITORY" => "FourierFlows/GeophysicalFlowsDocumentation") do
    example_filepath = joinpath(EXAMPLES_DIR, example)
    withenv("JULIA_DEBUG" => "Literate") do
      Literate.markdown(example_filepath, OUTPUT_DIR;
                        flavor = Literate.DocumenterFlavor(), execute = true)
      Literate.notebook(example_filepath, OUTPUT_DIR, execute = false)
      Literate.script(example_filepath, OUTPUT_DIR)
    end
  end
end

#####
##### Build and deploy docs
#####

format = Documenter.HTML(
  collapselevel = 2,
     prettyurls = get(ENV, "CI", nothing) == "true",
      canonical = "https://fourierflows.github.io/GeophysicalFlowsDocumentation/stable/"
)

makedocs(
 modules = [GeophysicalFlows],
 doctest = true,
   clean = true,
checkdocs = :all,
  format = format,
 authors = "Navid C. Constantinou, Gregory L. Wagner, and contributors",
sitename = "GeophysicalFlows.jl",
   pages = Any[
            "Home"    => "index.md",
            "Installation instructions" => "installation_instructions.md",
            "Aliasing" => "aliasing.md",
            "GPU" => "gpu.md",
            "Visualize output" => "visualize.md",
            "Examples" => [
              "TwoDNavierStokes" => Any[
                "literated/twodnavierstokes_decaying.md",
                "literated/twodnavierstokes_stochasticforcing.md",
                "literated/twodnavierstokes_stochasticforcing_budgets.md",
                ],
              "SingleLayerQG" => Any[
                "literated/singlelayerqg_betadecay.md",
                "literated/singlelayerqg_betaforced.md",
                "literated/singlelayerqg_decaying_topography.md",
                "literated/singlelayerqg_decaying_barotropic_equivalentbarotropic.md"
                ],
              "BarotropicQGQL" => Any[
                "literated/barotropicqgql_betaforced.md",
                ],
              "MultiLayerQG" => Any[
                "literated/multilayerqg_2layer.md"
                ],
              "SurfaceQG" => Any[
                "literated/surfaceqg_decaying.md"
                ]
            ],
            "Modules" => Any[
              "modules/twodnavierstokes.md",
              "modules/singlelayerqg.md",
              "modules/barotropicqgql.md",
              "modules/multilayerqg.md",
              "modules/surfaceqg.md"
            ],
            "Stochastic forcing" => "stochastic_forcing.md",
            "Contributor's guide" => "contributing.md",
            "Library" => Any[
              "lib/types.md",
              "lib/functions.md"
            ]
           ]
)

@info "Cleaning up temporary .jld2 and .nc files created by doctests..."

for file in vcat(glob("docs/*.jld2"), glob("docs/*.nc"))
    rm(file)
end

withenv("GITHUB_REPOSITORY" => "FourierFlows/GeophysicalFlowsDocumentation") do
  deploydocs(       repo = "github.com/FourierFlows/GeophysicalFlowsDocumentation.git",
                versions = ["stable" => "v^", "v#.#.#", "dev" => "dev"],
            push_preview = false,
               forcepush = true,
               devbranch = "main"
            )
end

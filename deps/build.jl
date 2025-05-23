using CondaPkg
include("config_init.jl")

@info "installing cif2cell"
CondaPkg.update()

if any(x -> !ispath(joinpath(@__DIR__, x)),
       ("wannier90flags.jl", "qeflags.jl", "abinitflags.jl", "elkflags.jl",
        "qe7.2flags.jl"))
    include("asset_init.jl")
end

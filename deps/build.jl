using CondaPkg
include("config_init.jl")

@info "installing cif2cell"
CondaPkg.add("cif2cell")

CondaPkg.withenv() do
  cif2cell = CondaPkg.which("cif2cell")
  run(`$cif2cell --version`)
end

if any(x->!ispath(joinpath(@__DIR__, x)), ("wannier90flags.jl", "qeflags.jl", "abinitflags.jl", "elkflags.jl", "qe7.2flags.jl"))
    include("asset_init.jl")
end

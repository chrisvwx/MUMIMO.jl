# --------------
# Initialization
# --------------
if Pkg.installed("LLLplus")==nothing
    Pkg.clone("git@github.com:christianpeel/LLLplus.jl.git")
end
Pkg.add("PyPlot")

include("../src/MUMIMO.jl")
using MUMIMO
using PyPlot


# --------------
# Basic test
# --------------
mimoUplink(100,4,[4],[4],[0],[12],12,[-5.0:5:20;])

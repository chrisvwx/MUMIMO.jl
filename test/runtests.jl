# --------------
# Initialization
# --------------
Pkg.add("LLLplus")
Pkg.add("PyPlot")
if VERSION<=v"0.6.3"
    include("../src/MUMIMO.jl")
else
    using Pkg
end

using MUMIMO
using PyPlot


# --------------
# Basic test
# --------------
mimoUplink(100,4,[4],[4],[0],[12],12,[-5.0:5:20;])

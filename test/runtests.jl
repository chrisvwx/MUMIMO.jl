# --------------
# Initialization
# --------------
if VERSION<=v"0.6.3"
    include("../src/MUMIMO.jl")
else
    using Pkg
end
Pkg.add("LLLplus")
Pkg.add("PyPlot")

using MUMIMO
using PyPlot


# --------------
# Basic test
# --------------
#mimoUplink(1,4,[2],[2],[0],[12],1,[30.0;])
mimoUplink(100,4,[4],[4],[0],[12],12,[-5.0:5:20;])

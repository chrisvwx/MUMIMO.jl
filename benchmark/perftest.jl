Pkg.add("LLLplus")

#include("../src/MUMIMO.jl")
using MUMIMO
using Plots

mimoUplink(10000,4,[4],[4],[0],[12],12,[-5.0:5:20;])
savefig("benchmark/perfVsSNRqpsk4ant.png") # run from root directory
mimoUplink(10000,4,[8],[8],[0],[12],12,[-5.0:5:20;])
savefig("benchmark/perfVsSNRqpsk8ant.png") # run from root directory


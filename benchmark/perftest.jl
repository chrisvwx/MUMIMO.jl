if Pkg.installed("LLLplus")==Nothing
    Pkg.clone("git@github.com:christianpeel/LLLplus.jl.git")
end
Pkg.add("PyPlot")

include("../src/MUMIMO.jl")
using MUMIMO
using PyPlot

mimoUplink(1000,4,[4],[4],[0],[12],12,[-5.0:5:20;])
savefig("benchmark/perfVsSNRqpsk4ant.png") # run from root directory
mimoUplink(1000,4,[8],[8],[0],[12],12,[-5.0:5:20;])
savefig("benchmark/perfVsSNRqpsk8ant.png") # run from root directory
mimoUplink(1000,4,[4],[3],[1],[12],12,[10.0;],[-15.0:5:15])
savefig("benchmark/perfVsCIRqpsk4ant.png") # run from root directory
mimoUplink(1000,4,[2:8],[2],[0],[12],12,[10.0])
savefig("benchmark/perfVsMqpsk2user.png") # run from root directory


module MUMIMO

using LLLplus
using PyPlot

export
    mimoUplink,
    scodes,
    mimo_slice,
    siso_demod

include("mimoUplink.jl")
include("scodes.jl")

end # LLLplus module

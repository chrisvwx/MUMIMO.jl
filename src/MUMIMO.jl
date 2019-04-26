module MUMIMO

using LinearAlgebra
using LLLplus
using Plots
using Printf
using Statistics

export
    mimoUplink,
    scodes,
    mimo_slice,
    siso_demod

include("mimoUplink.jl")
include("scodes.jl")

end # LLLplus module

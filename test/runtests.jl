# --------------
# Initialization
# --------------
Pkg.add("LLLplus")
using Pkg
using MUMIMO

# --------------
# Basic test
# --------------
mimoUplink(100,4,[4],[4],[0],[12],12,[-5.0:5:20;])

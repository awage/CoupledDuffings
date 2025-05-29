using DrWatson
@quickactivate 
using CairoMakie
using LaTeXStrings
using Attractors

include(srcdir("duffing_mapper.jl"))

# N = 5; c = 0.1; k1 = 1.; k3 = 0.1; F = 0.4; kc = 0.05; ω = 1.2457
N = 3; c = 0.1; k1 = 1.0; k3 = 0.1; F = 0.4; kc = 0.05; ω = 1.2457
# prange = range(0.05, 0.5, length = 50)
prange = range(1.3, 1.8, length = 10)
# prange = range(0.05, 0.5, step = 0.01)
dps = duffing_parameters(N, c, k1, k3, F, kc, ω) 

# compute globl continuation of attractors
yg = range(-15,15, length=10)
grid = ntuple(x -> yg, N*2)
sampler, _ = statespace_sampler(grid) 
mapper = get_mapper(dps)

ascm = AttractorSeedContinueMatch(mapper)
pidx = :ω

fractions_cont, attractors_cont = global_continuation(
	ascm, prange, pidx, sampler; samples_per_parameter = 2000
)

fig = plot_basins_attractors_curves(
	fractions_cont, attractors_cont, A -> minimum(A[:, 1]), prange,
)

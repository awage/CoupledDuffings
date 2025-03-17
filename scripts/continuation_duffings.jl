using DrWatson
@quickactivate 
using CairoMakie
using LaTeXStrings
using Attractors

include(srcdir("duffing_mapper.jl"))

N = 5; c = 0.1; k1 = 1.; k3 = 0.1; F = 0.4; kc = 0.05; ω = 1.2457
dps = duffing_parameter(N, c, k1, k3, F, kc, ω) 

# compute globl continuation of attractors
yg = range(-5,5, length=10)
grid = ntuple(x -> yg, N*2)
sampler, _ = statespace_sampler(grid) 
mapper = get_mapper(dps)

ascm = AttractorSeedContinueMatch(mapper)
prange = range(0.01, 0.1, length = 10)
pidx = :kc

fractions_cont, attractors_cont = global_continuation(
	ascm, prange, pidx, sampler; samples_per_parameter = 1000
)

fig = plot_basins_attractors_curves(
	fractions_cont, attractors_cont, A -> minimum(A[:, 1]), prange,
)

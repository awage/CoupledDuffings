using DrWatson
@quickactivate 
using CairoMakie
using LaTeXStrings
using Attractors

include(srcdir("duffing_mapper.jl"))


function compute_k_sweep(params::Dict)
    @unpack  Np, Nsamples = params

    N = 3; c = 0.1; k1 = 1.0; k3 = 0.1; F = 0.4; kc = 0.05; ω = 1.2457
    prange = range(0.05, 0.5, length = Np)
    dps = duffing_parameters(N, c, k1, k3, F, kc, ω) 

    # compute globl continuation of attractors
    yg = range(-5,5, length=10)
    grid = ntuple(x -> yg, N*2)
    sampler, _ = statespace_sampler(grid) 
    mapper = get_mapper(dps)

    matcher = MatchBySSSetDistance(; distance = Hausdorff())
    ascm = AttractorSeedContinueMatch(mapper, matcher)
    pidx = :kc

    fractions_cont, attractors_cont = global_continuation(
        ascm, prange, pidx, sampler; samples_per_parameter = Nsamples
    )

    return @strdict(fractions_cont,  attractors_cont, prange, dps)
end


Np = 5; Nsamples = 500
params = @strdict Np Nsamples

dat, _ = produce_or_load(compute_k_sweep, params; prefix = "k_sweep", force = false)

@unpack fractions_cont, attractors_cont, prange = dat

fig = plot_basins_attractors_curves(
	fractions_cont, attractors_cont, A -> minimum(A[:, 1]), prange,
)


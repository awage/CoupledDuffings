using DrWatson
@quickactivate 
using CairoMakie
using Statistics
using LaTeXStrings
using Attractors

include(srcdir("duffing_mapper.jl"))

function compute_freq_sweep(params::Dict)
    @unpack  Np, Nsamples, k3, kc = params
    N = 5; c = 0.1; k1 = 1.0;  F = 0.4;  ω = 1.2457
    # prange = range(0.05, 0.5, length = 50)
    prange = range(1., 2., length = Np)
    # prange = range(0.05, 0.5, step = 0.01)
    dps = duffing_parameters(N, c, k1, k3, F, kc, ω) 

    # compute globl continuation of attractors
    yg = range(-5,5, length=10)
    grid = ntuple(x -> yg, N*2)
    sampler = RandomICsSampler(Nsamples, grid) 
    mapper = get_mapper(dps)

    matcher = MatchBySSSetDistance(; distance = Hausdorff())
    ascm = AttractorSeedContinueMatch(mapper, matcher)
    pidx = :ω
    pcurve = [Dict(pidx => p) for p in prange]
    fractions_cont, attractors_cont = global_continuation(ascm, pcurve, sampler)

    return @strdict(fractions_cont,  attractors_cont, prange, dps)
end


k3 = 0.5; 
kc = 0.3;
Np = 15; Nsamples = 1000; N = 5
params = @strdict Np Nsamples kc k3

dat, _ = produce_or_load(compute_freq_sweep, params; prefix = "freq_sweep", force = false)
@unpack fractions_cont, attractors_cont, prange = dat

# matcher = MatchBySSSetDistance(; distance = Hausdorff())
# rmaps = match_sequentially!(attractors_cont, matcher)
# match_sequentially!(fractions_cont, rmaps)

fig = plot_basins_attractors_curves(fractions_cont, attractors_cont, A -> maximum(mean(Matrix(A),dims = 2)), prange)
# fig = plot_basins_attractors_curves(fractions_cont, attractors_cont, A -> sum(mean(Matrix(A),dims = 2)), prange)
# fig = plot_basins_attractors_curves(fractions_cont, attractors_cont, A -> maximum(std(Matrix(A),dims = 2)), prange)

# @save "freq_sweep.jld2" fractions_cont attractors_cont fig

using DrWatson
@quickactivate 
using CairoMakie
using LaTeXStrings
using Attractors
using OrdinaryDiffEq:Vern9
using ProgressMeter
using LinearAlgebra 

include(srcdir("duffing_mapper.jl"))

function compute_basins_random(dps::DuffingParameters)
    #random ics for now
    N = dps.N
    mapper = get_mapper(dps)
    Nics = 100
    bsn = @showprogress [ mapper(4*(rand(N*2).- 0.5)) for k in 1:Nics]
    att = extract_attractors(mapper)
    return @strdict(bsn,  att)
end

N = 5; c = 0.1; k1 = 1.; k3 = 0.1; F = 0.4; kc = 0.05; ω = 1.2457
dps = duffing_parameters(N, c, k1, k3, F, kc, ω) 

# compute basins



dat, _ = produce_or_load(compute_basins_random, dps; prefix = "coupled_duffings", force = false)

@unpack bsn, att = dat

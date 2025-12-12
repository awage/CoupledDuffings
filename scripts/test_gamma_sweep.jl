using DrWatson
@quickactivate 
using CairoMakie
using Statistics
using LaTeXStrings
using Attractors

include(srcdir("duffing_mapper.jl"))

function compute_gamma_sweep(params::Dict)
    @unpack  N, Np, Nsamples, k3, ω = params
     c = 0.1; k1 = 1.0;  F = 0.4; kc =0.1 #ω = 1.2457
    prange = 10 .^range(-1,0,length = Np)
    # \gamma = kc/k3
    kc_range = k3.*prange
    
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
        ascm, kc_range, pidx, sampler; samples_per_parameter = Nsamples
    )

    return @strdict(fractions_cont,  attractors_cont, prange, dps)
end

function get_number_attractors(k3, ω)
    kc = 0.3;
    Np = 40; Nsamples = 1000; 

    Natt = Vector{Vector{Float64}}()
    γrange = []
    Threads.@threads for N in [2 3 4 5]

        params = @strdict Np Nsamples k3 ω N

        dat, _ = produce_or_load(compute_gamma_sweep, params; prefix = "inv_gamma_sweep", force = false)
        @unpack fractions_cont, attractors_cont, prange = dat
        γrange = prange
        v = zeros(Np)
        for (k,a) in enumerate(attractors_cont)
            v[k] = length(a)
        end

        push!(Natt, v)
    end

    return Natt, γrange
end

k3 = 0.1; 
ω = 1.2457
Natt, γrange = get_number_attractors(k3, ω)
f = Figure()
ax = Axis(f[1,1], xscale = log10, yscale = log2,  yticks = ([2, 4, 8, 16, 32, 1024], ["2", "4", "8", "16", "32", "1024"]), 
      xticks = ([1, 2, 5, 10], ["1", "2", "5", "10"]), 
      yticklabelsize = 20,  ylabelsize = 25, xticklabelsize = 20,  xlabelsize = 25,  xminorticksvisible = true, xminorgridvisible = true,
      titlesize = 20,   xminorticks = IntervalsBetween(5), xlabel = L"\gamma", ylabel = L"N_{att}")
for v in Natt 
    lines!(ax, γrange, v)
end

save("Natt_k3=0.1.pdf", f)

k3 = 0.5; 
ω = 1.5035
Natt, γrange = get_number_attractors(k3, ω)
f = Figure()
ax = Axis(f[1,1], xscale = log10, yscale = log2,  yticks = ([2, 4, 8, 16, 32, 1024], ["2", "4", "8", "16", "32", "1024"]), 
      # xticks = ([1, 2, 5, 10], ["1", "2", "5", "10"]), 
      yticklabelsize = 20,  ylabelsize = 25, xticklabelsize = 20,  xlabelsize = 25,  xminorticksvisible = true, xminorgridvisible = true,
      titlesize = 20,   xminorticks = IntervalsBetween(5), xlabel = L"\gamma", ylabel = L"N_{att}")
for v in Natt 
    lines!(ax, γrange, v)
end

save("Natt_k3=0.5.pdf", f)

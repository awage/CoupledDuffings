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
    prange = 10 .^range(-1,log10(1.1),length = Np)
    # \gamma = kc/k3
    kc_range = k3.*prange
    
    dps = duffing_parameters(N, c, k1, k3, F, kc, ω) 

    # compute globl continuation of attractors
    yg = range(-5,5, length=10)
    grid = ntuple(x -> yg, N*2)
    # Uniform random ics in the grid box. RandomICsSampler(Nsamples, grid) is not used because
    # statespace_sampler indexes per-thread buffers by threadid() and fails inside Threads.@threads.
    sampler = RandomICsSampler(() -> first(yg) .+ (last(yg) - first(yg)) .* rand(N*2), Nsamples)
    mapper = get_mapper(dps)

    matcher = MatchBySSSetDistance(; distance = Hausdorff())
    ascm = AttractorSeedContinueMatch(mapper, matcher)
    pcurve = [Dict(:kc => kc) for kc in kc_range]
    fractions_cont, attractors_cont = global_continuation(ascm, pcurve, sampler)

    return @strdict(fractions_cont,  attractors_cont, prange, dps)
end

# Attractor dicts saved with an older StateSpaceSet are read back by JLD2 as SerializedDict
as_dict(a) = a isa AbstractDict ? a : Dict(a.kvvec)

# Amplitude of each oscillator estimated from the stroboscopic points: A_i = sqrt(x_i^2 + (v_i/ω)^2)
function oscillator_amplitudes(A, N, ω)
    P = A isa StateSpaceSet ? vec(A) : A.data # A.data: StateSpaceSet reconstructed by JLD2
    return [maximum(sqrt(p[i]^2 + (p[N+i]/ω)^2) for p in P) for i in 1:N]
end

# Uniform response (L^N or H^N): all oscillators share the same amplitude level.
# Low and high branches are well separated (≈0.8 vs ≈2.8), so a relative spread criterion suffices.
function is_uniform(A, N, ω; δ = 0.25)
    amp = oscillator_amplitudes(A, N, ω)
    return (maximum(amp) - minimum(amp)) < δ*maximum(amp)
end

# P_U(γ) = sum of basin fractions of the uniform attractors {L^N, H^N}
function uniform_fraction(attractors_cont, fractions_cont, N, ω)
    PU = zeros(length(attractors_cont))
    for (k, a) in enumerate(attractors_cont)
        PU[k] = sum(get(fractions_cont[k], key, 0.0) for (key, A) in as_dict(a) if is_uniform(A, N, ω); init = 0.0)
    end
    return PU
end

function get_number_attractors(k3, ω)
    kc = 0.3;
    Np = 30; Nsamples = 10000; 

    Nrange = [2, 3, 4, 5, 10]
    Natt = Vector{Vector{Float64}}(undef, length(Nrange))
    PU = Vector{Vector{Float64}}(undef, length(Nrange))
    γrange = []
    Threads.@threads for j in eachindex(Nrange)
        N = Nrange[j]

        params = @strdict Np Nsamples k3 ω N

        dat, _ = produce_or_load(compute_gamma_sweep, params; prefix = "inv_gamma_sweep", force = false)
        @unpack fractions_cont, attractors_cont, prange = dat
        γrange = prange
        v = zeros(Np)
        for (k,a) in enumerate(attractors_cont)
            v[k] = length(as_dict(a))
        end

        Natt[j] = v
        PU[j] = uniform_fraction(attractors_cont, fractions_cont, N, ω)
    end

    return Natt, PU, γrange
end

function plot_Natt_PU(k3, ω, fname)
    Natt, PU, γrange = get_number_attractors(k3, ω)
    xargs = (xscale = log10, xticks = ([0.1, 1], ["0.1", "1"]), 
          yticklabelsize = 20,  ylabelsize = 25, xticklabelsize = 20,  xlabelsize = 25,  xminorticksvisible = true, xminorgridvisible = true,
          titlesize = 20,   xminorticks = IntervalsBetween(5))
    lab =["N = 2", "N = 3", "N = 4", "N = 5", "N = 10"] 

    f = Figure(size = (600, 700))
    ax = Axis(f[1,1]; yscale = log2,  yticks = ([2, 4, 8, 16, 32, 1024], ["2", "4", "8", "16", "32", "1024"]), 
          ylabel = L"N_{att}", xargs...)
    for (k,v) in enumerate(Natt) 
        lines!(ax, γrange, v; label = lab[k])
    end
    axislegend(ax)
    hidexdecorations!(ax, grid = false, minorgrid = false, ticks = false, minorticks = false)

    ax2 = Axis(f[2,1]; xlabel = L"\gamma", ylabel = L"P_{\mathrm{U}}", xargs...)
    for (k,v) in enumerate(PU) 
        lines!(ax2, γrange, v; label = lab[k])
    end
    ylims!(ax2, 0, 1.05)
    linkxaxes!(ax, ax2)

    save(fname, f)
end

plot_Natt_PU(0.1, 1.2457, "Natt_k3=0.1.pdf")
plot_Natt_PU(0.5, 1.5035, "Natt_k3=0.5.pdf")

using DrWatson
@quickactivate 
using CairoMakie
using ProgressMeter
include(srcdir("duffing_mapper.jl"))


function compute_basins(mapper, N)
    #random ics for now
    Nics = 10000
    xg = yg = range(-5,5, length = 100) 
    bsn = @showprogress [ mapper([x; 0.; y; 0.; zeros(2*(N-2))]) for x in xg, y in yg]
    att = extract_attractors(mapper)
    return bsn,  att
end

N = 9; c = 0.1; k1 = 1.0; k3 = 0.5; F = 0.4; kc = 0.5; ω = 1.5035
dps = duffing_parameters(N, c, k1, k3, F, kc, ω) 


# compute basins
mapper = get_mapper(dps)
bsn, att = compute_basins(mapper, N)

xg = yg = range(-5,5, length = 100) 
heatmap(xg, yg, bsn)

using DrWatson
@quickactivate 
using CairoMakie
using LaTeXStrings
using Attractors
using OrdinaryDiffEq:Vern9
using ProgressMeter
using LinearAlgebra 


function coupled_duffings!(du, u, p, t)
    N = p[1];  W = p[2]
    c, k1, k3, F, kc, ω = p[3:8]

    uu = u[1:N]
    duu = u[N+1:2*N]

    cpld = kc*W*uu # Coupling term

    du[1:N] = duu
    @. du[N+1:end] .= -c*duu + -k1*uu - k3*uu^3 - cpld +  F*sin(ω*t)
    return nothing
end

# coupling with first neighbors
function coupling_matrix(N)
    W = diagm(1 => -ones(N-1)) + diagm(-1 => -ones(N-1))
    rowsum = W*ones(N)
    W = W + diagm(0 => -rowsum) 
end

function compute_basins(di::Dict)
    @unpack  N, c, k1, k3, F, kc, ω = di
    W = coupling_matrix(N)
    @show W
    diffeq = (alg = Vern9(), reltol = 1e-10, maxiters = 1e8)
    ds = CoupledODEs(coupled_duffings!, rand(N*2), (N, W, c, k1, k3, F, kc, ω); diffeq)
    smap = StroboscopicMap(ds, 2*pi/ω) # Stroboscopic map definition
    yg =  collect(range(-10, 10; length = 50001))
    grid = ntuple(x -> yg, N*2)
    mapper = AttractorsViaRecurrences(smap, grid; 
                    consecutive_basin_steps = 100, 
                    consecutive_recurrences = 1000,
                    attractor_locate_steps = 1000)
    #random ics for now
    Nics = 10000
    xg = yg = range(-5,5, length = 100) 
    # bsn = @showprogress [ mapper(4*(rand(N*2).- 0.5)) for k in 1:Nics]
    bsn = @showprogress [ mapper([x; 0.; y; 0.; zeros(2*(N-2))]) for x in xg, y in yg]

    att = extract_attractors(mapper)
    return @strdict(bsn,  att)
end

function get_trajectory(u0, T, di::Dict)
    @unpack  N, c, k1, k3, F, kc, ω = di
    W = coupling_matrix(N)
    diffeq = (alg = Vern9(), reltol = 1e-6, maxiters = 1e8)
    ds = CoupledODEs(coupled_duffings!, rand(N*2), (N, W, c, k1, k3, F, kc, ω); diffeq)
    smap = StroboscopicMap(ds, 2*pi/ω) # Stroboscopic map definition
    y,t = trajectory(smap, T, u0) 
    return  y,t 
end

# N = 5; c = 0.1; k1 = 1.; k3 = 0.1; F = 0.4; kc = 0.05; ω = 1.2457
N = 9; c = 0.1; k1 = 1; k3 = 0.5; F = 0.4; kc = 0.5; ω = 1.5035

params = @strdict N c k1 k3 F kc ω

# get a sample trajectory: 
# y,t = get_trajectory(rand(2*N), 100, params)
# scatter(t, y[:,1])


# compute basins
dd = compute_basins(params)
@unpack bsn, att = dd

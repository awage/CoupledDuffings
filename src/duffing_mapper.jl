using Attractors
# using OrdinaryDiffEq:Vern9
using OrdinaryDiffEqVerner
using LinearAlgebra 


mutable struct DuffingParameters{N}
    N::Int
    W::Array
    c::N
    k1::N
    k3::N
    F::N
    kc::N
    ω::N
end

function duffing_parameters(N, c, k1, k3, F, kc, ω)
    W = coupling_matrix(N)
    return DuffingParameters(N, W, c, k1, k3, F, kc, ω)
end

function coupled_duffings!(du, u, p, t)
    (;N, W, c, k1, k3, F, kc, ω) = p

    uu = u[1:N]
    duu = u[N+1:2*N]

    cpld = kc*W*uu # Coupling term

    du[1:N] = duu
    @. du[N+1:end] .= -c*duu + -k1*uu - k3*uu^3 - cpld +  F*cos(ω*t)
    return nothing
end

# coupling with first neighbors
function coupling_matrix(N)
    W = diagm(1 => -ones(N-1)) + diagm(-1 => -ones(N-1))
    rowsum = W*ones(N)
    W = W + diagm(0 => -rowsum) 
end

function get_mapper(dps::DuffingParameters)
    (;N, W, c, k1, k3, F, kc, ω) = dps
    diffeq = (alg = Vern9(), reltol = 1e-8, maxiters = 1e6)
    ds = CoupledODEs(coupled_duffings!, rand(N*2), dps; diffeq)
    smap = StroboscopicMap(ds, 2*pi/ω) # Stroboscopic map definition
    yg =  collect(range(-20, 20; length = 5001))
    grid = ntuple(x -> yg, N*2)
    mapper = AttractorsViaRecurrences(smap, grid; 
                    consecutive_basin_steps = 100, 
                    consecutive_recurrences = 2000,
                    attractor_locate_steps = 2000)
    return mapper
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

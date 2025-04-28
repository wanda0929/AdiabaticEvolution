using Test
using Plots
using LinearAlgebra
using Yao, Yao.EasyBuild
using AdiabaticEvolution

#using DifferentialEquations

@testset "ARP_pulse" begin
    Ω = 2.0 * 2π
    α = 1.0 * 2π
    τ = 1.0
    
    t_total = 16.0
    times = -8:0.01:8
    
    # Calculate values for both functions
    Omega_0 = Rabi_ARP(Ω, τ, α)
    Delta_0 = detuning_ARP(α)
    Omegavec = [Omega_0(t) for t in times]
    Deltavec = [Delta_0(t) for t in times]
    
    # Create the plot
    p = Plots.plot(times, [real(Omegavec)/2/pi, real(Deltavec)/2/pi], 
        label=["Rabi Frequency Ω(t)" "Detuning Δ(t)"],
        title="ARP Pulse Parameters",
        xlabel="Time",
        ylabel="Value",
        linewidth=2)
    display(p)
    @test p isa Plots.Plot
end

@testset "ARP_hamiltonian" begin
    Ω = 2.0
    α = 1.0
    τ = 1.0
    U = 100.0
    t_total = 16.0
    times = 0:0.01:t_total
    h = ARP_hamiltonian_2atoms(Ω, τ, α, U)
    @test h isa Function
    #Matrix(h(1.0))
end

@testset "two-atom phase accumulation" begin
    Ω = 2.0 * 2π
    α = 1.0 * 2π
    τ = 1.0
    U = 20.0 * 2π
    t_total = 16.0
    reg = zero_state(2)
    h = ARP_hamiltonian_2atoms(Ω, τ, α, U)
    reg, times, phases = evolve_ARP(reg, h, t_total)
    p = Plots.plot(times, phases, 
    label=["Phase"],
    title="two-atom phase accumulation",
    xlabel="Time",
    ylabel="Phase",
    linewidth=2)
    display(p)
    @test times isa Vector{Float64}
    @test phases isa Vector{Float64}
end
@testset "single-atom phase accumulation" begin 
    Ω = 2.0 * 2π
    α = 1.0 * 2π
    τ = 1.0
    t_total = 16.0
    reg = zero_state(1)
    h = ARP_hamiltonian_singleatom(Ω, τ, α)
    reg, times, phases = evolve_ARP(reg, h, t_total)
    p = Plots.plot(times, phases, 
    label=["Phase"],
    title="Single-atom phase accumulation",
    xlabel="Time",
    ylabel="Phase",
    linewidth=2)
    display(p)
    @test times isa Vector{Float64}
    @test phases isa Vector{Float64}
end


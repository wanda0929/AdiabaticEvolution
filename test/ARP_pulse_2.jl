using Test
using Plots
using LinearAlgebra
using Yao, Yao.EasyBuild
using AdiabaticEvolution
import Yao: expect
#using DifferentialEquations

@testset "ARP_pulse" begin
    Ω = 2.0 * 2π
    α = 1.0 * 2π
    τ = 1.0
    
    t_total = 8.0
    times = -4:0.01:4
    
    # Calculate values for both functions
    Omega_0 = Rabi_ARP_2(Ω, τ, α)
    Delta_0 = detuning_ARP_2(α)
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
    Ω = 2.0 * 2π
    α = 1.0 * 2π
    τ = 1.0
    U = 1000.0 * 2π
    t_total = 16.0
    times = 0:0.01:t_total
    h = ARP_hamiltonian_2atoms(Ω, τ, α, U)
    @test h isa Function
    #Matrix(h(1.0))
end

@testset "two-atom phase accumulation" begin
    Ω1 = 2.0 * 2π
    Ω2 = 0.5 * 2π
    α1 = 1.0 * 2π
    α2 = 0.5 * 2π
    τ = 1.0
    U = 2000.0 * 2π
    t_total = 8.0
    reg1 = zero_state(2)
    reg2 = zero_state(2)
    h1 = ARP_hamiltonian_2atoms(Ω1, τ, α1, U)
    h2 = ARP_hamiltonian_2atoms(Ω2, τ, α2, U)
    reg1, times, phases1, populations1 = evolve_ARP_2(reg1, h1, t_total)
    reg2, times, phases2, populations2 = evolve_ARP_2(reg2, h2, t_total)
    p = Plots.plot(times, [populations1, populations2], 
    label=["Population Ω1/2π=2, α1/2π=1" "Population Ω2/2π=0.5, α2/2π=0.5"],
    title="two-atom population",
    xlabel="Time",
    ylabel="Population",
    linewidth=2)

    #p1 = Plots.plot(times, phases, 
    #label=["Phase"],
    #title="two-atom phase accumulation",
    #xlabel="Time",
    #ylabel="Phase",
    #linewidth=2)
    
    display(p)
    #display(p1)
    @test times isa Vector{Float64}
    @test phases isa Vector{Float64}
end



@testset "single-atom evolution" begin
    Ω1 = 2.0 * 2π
    Ω2 = 0.5 * 2π
    α1 = 1.0 * 2π
    α2 = 0.5 * 2π
    τ = 1.0
    t_total = 8.0
    reg1 = zero_state(1)
    reg2 = zero_state(1)
    h1 = ARP_hamiltonian_singleatom_2(Ω1, τ, α1)
    h2 = ARP_hamiltonian_singleatom_2(Ω2, τ, α2) 
    reg1, times, phases1, populations1 = evolve_ARP_2(reg1, h1, t_total)
    reg2, times, phases2, populations2 = evolve_ARP_2(reg2, h2, t_total)
    
    # 绘制相位
    p1 = Plots.plot(times, phases1, 
        label="Phase Ω=2π",
        title="Single-atom Phase Evolution",
        xlabel="Time",
        ylabel="Phase",
        linewidth=2)
    
    p4 = Plots.plot!(times, phases2,
        label="Phase Ω=0.5π",
        linewidth=2)
    
    # 绘制布居数
    p2 = Plots.plot(times, [populations1, populations2], 
        label=["Population Ω/2π=2, α/2π=1" "Population Ω/2π=0.5, α/2π=0.5"],
        title="Single-atom Population Evolution",
        xlabel="Time",
        ylabel="Population",
        linewidth=2)
    

    display(p2)
    
    @test times isa Vector{Float64}
    @test phases1 isa Vector{Float64}
    @test populations1 isa Vector{Float64}
end


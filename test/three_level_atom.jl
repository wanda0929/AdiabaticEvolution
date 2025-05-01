using Test
using Plots
using LinearAlgebra
using QuantumOptics
using Revise
using AdiabaticEvolution
using OrdinaryDiffEq

@testset "STIRAP pulse" begin
    Ω1 = 10.0 * 2π
    Ω2 = 10.0 * 2π
    τ = 1.0
    Dt = 2.0
    pulse_s = STIRAP_pulse_s(Ω1, τ)
    pulse_c = STIRAP_pulse_c(Ω2, τ, Dt)
    times = range(0, 12, length=1000)
    Omegavec = [pulse_s(t) for t in times]
    Deltavec = [pulse_c(t) for t in times]
    
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



@testset "Three-level atom Hamiltonian" begin
    # Test parameters
    Ω1 = 10.0 * 2π
    Ω2 = 10.0 * 2π
    τ = 1.0
    Dt = 2.0
    Ω1_func = STIRAP_pulse_c(Ω1, τ, Dt)
    Ω2_func = STIRAP_pulse_s(Ω2, τ)
    Δ1 = 0.0
    Δ2 = 0.0
    
    # Test Hamiltonian
    H_func = three_level_hamiltonian_time_dependent(
        Ω1_func,  # Constant Ω1
        Ω2_func,  # Constant Ω2
        Δ1,
        Δ2
    )
    
    # Test at specific time
    H = H_func(1.0, create_initial_state(1))
    @test isa(H, Operator)
    @test length(H.basis_l) == 3
    @test length(H.basis_r) == 3
    
    # Test evolution
    ψ0 = create_initial_state(1)  # Start in state |1⟩
    t_total = 10.0
    final_ψ, times, populations = evolve_three_level(ψ0, H_func, t_total)
    
    # 绘制布居数演化
    p = plot(times, [pop[1] for pop in populations], label="|1⟩")
    plot!(p, times, [pop[2] for pop in populations], label="|2⟩")
    plot!(p, times, [pop[3] for pop in populations], label="|3⟩")
    xlabel!("Time")
    ylabel!("Population")
    title!("STIRAP evolution")
display(p)
end 


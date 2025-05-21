using Test
using Plots
using LinearAlgebra
using QuantumOptics
using Revise
using AdiabaticEvolution
using OrdinaryDiffEq
using Yao, Yao.EasyBuild

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

@testset "Three-level atom dynamics" begin
    # 设置参数
    Ω = 1.0  # Rabi频率
    τ = 1.0  # 脉冲宽度
    Dt = 2.0  # 脉冲延迟
    t_total = 10.0  # 总演化时间
    
    # 创建STIRAP脉冲
    Ω1_func, Ω2_func, Δ1_func, Δ2_func = create_stirap_pulses(Ω, τ, Dt)
    
    # 创建哈密顿量
    H = create_time_dependent_hamiltonian(Ω1_func, Ω2_func, Δ1_func, Δ2_func)
    
    # 创建初始态
    reg = create_initial_state()
    
    # 演化系统
    reg, times, populations = evolve_three_level(reg, H, t_total)
    
    # 绘制布居数演化
    p = plot(times, hcat(populations...), 
        label=["|g⟩" "|e⟩" "|r⟩"],
        title="Three-level atom population evolution",
        xlabel="Time",
        ylabel="Population",
        linewidth=2)
    display(p)
    
    # 测试
    @test length(times) == length(populations)
    @test length(populations[1]) == 3  # 三个能级的布居数
    @test sum(populations[end]) ≈ 1.0  # 布居数守恒
end


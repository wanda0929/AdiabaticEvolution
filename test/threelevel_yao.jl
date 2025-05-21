using Test
using Plots
using LinearAlgebra
using Yao, Yao.EasyBuild
using AdiabaticEvolution

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
        title="STIRAP Pulse Parameters",
        xlabel="Time",
        ylabel="Value",
        linewidth=2)
    display(p)
    @test p isa Plots.Plot
end

@testset "Three-level atom dynamics" begin
    # 设置参数
    Ω1 = 10.0 * 2π
    Ω2 = 10.0 * 2π
    τ = 1.0
    Dt = 2.0
    Δ1 = 0.0
    Δ2 = 0.0
    t_total = 12.0
    
    # 创建STIRAP脉冲
    pulse_s = STIRAP_pulse_s(Ω1, τ)
    pulse_c = STIRAP_pulse_c(Ω2, τ, Dt)
    
    # 创建哈密顿量
    H = three_level_hamiltonian_time_dependent(pulse_s, pulse_c, Δ1, Δ2)
    
    # 创建初始态
    reg = create_initial_state()
    
    # 演化系统
    reg, times, populations = evolve_three_level(reg, H, t_total)
    
    # 绘制布居数演化
    p = plot(times, hcat(populations...), 
        label=["|1⟩" "|2⟩" "|3⟩"],
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

@testset "Two three-level qubits dynamics" begin
    # 设置参数
    Ω1 = 10.0 * 2π  # 第一个量子比特的0-1跃迁
    Ω2 = 10.0 * 2π  # 第一个量子比特的1-2跃迁
    Ω3 = 10.0 * 2π  # 第二个量子比特的0-1跃迁
    Ω4 = 10.0 * 2π  # 第二个量子比特的1-2跃迁
    τ = 1.0
    Dt = 2.0
    Δ1 = 0.0  # 第一个量子比特的|1⟩能级失谐
    Δ2 = 0.0  # 第一个量子比特的|2⟩能级失谐
    Δ3 = 0.0  # 第二个量子比特的|1⟩能级失谐
    Δ4 = 0.0  # 第二个量子比特的|2⟩能级失谐
    t_total = 12.0
    
    # 创建STIRAP脉冲
    pulse1 = STIRAP_pulse_s(Ω1, τ)  # 第一个量子比特的0-1跃迁
    pulse2 = STIRAP_pulse_s(Ω2, τ)  # 第一个量子比特的1-2跃迁
    pulse3 = STIRAP_pulse_s(Ω3, τ)  # 第二个量子比特的0-1跃迁
    pulse4 = STIRAP_pulse_s(Ω4, τ)  # 第二个量子比特的1-2跃迁
    
    # 创建哈密顿量
    H = three_level_hamiltonian_time_dependent(pulse1, pulse2, pulse3, pulse4, Δ1, Δ2, Δ3, Δ4)
    
    # 创建初始态
    reg = create_initial_state()
    
    # 演化系统
    reg, times, populations = evolve_three_level(reg, H, t_total)
    
    # 绘制布居数演化
    p = plot(times, hcat(populations...), 
        label=["Q1|0⟩" "Q1|1⟩" "Q1|2⟩" "Q2|0⟩" "Q2|1⟩" "Q2|2⟩"],
        title="Two three-level qubits population evolution",
        xlabel="Time",
        ylabel="Population",
        linewidth=2)
    display(p)
    
    # 测试
    @test length(times) == length(populations)
    @test length(populations[1]) == 6  # 两个量子比特，每个有三个能级
    @test sum(populations[end]) ≈ 1.0  # 布居数守恒
end

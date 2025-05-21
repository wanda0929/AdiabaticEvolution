using Yao, Yao.EasyBuild
using LinearAlgebra

# 导入张量积运算符
import Yao: ⊗

# 定义STIRAP脉冲
function STIRAP_pulse_s(Ω::Float64, τ::Float64)
    function O(t::Float64)
        return Ω * (exp(-(t-5)^2/(2*τ^2)))
    end
    return O
end

function STIRAP_pulse_c(Ω::Float64, τ::Float64, Dt::Float64)
    function O(t::Float64)
        return Ω * (exp(-(t-5-Dt)^2/(2*τ^2)))
    end
    return O
end

# 定义三能级量子比特的基态和激发态投影算符
function create_projectors()
    # 第一个量子比特的投影算符
    P1_1 = put(2, 1=>I2) * put(1, 1=>I2-Z)  # |0⟩⟨0|
    P1_2 = put(2, 1=>I2) * put(1, 1=>I2+Z)  # |1⟩⟨1|
    P1_3 = put(2, 1=>I2-Z) * put(1, 1=>I2)  # |2⟩⟨2|
    
    # 第二个量子比特的投影算符
    P2_1 = put(2, 2=>I2) * put(1, 2=>I2-Z)  # |0⟩⟨0|
    P2_2 = put(2, 2=>I2) * put(1, 2=>I2+Z)  # |1⟩⟨1|
    P2_3 = put(2, 2=>I2-Z) * put(1, 2=>I2)  # |2⟩⟨2|
    
    return P1_1, P1_2, P1_3, P2_1, P2_2, P2_3
end

# 定义三能级量子比特的跃迁算符
function create_transition_operators()
    # 第一个量子比特的跃迁算符
    S1_01 = put(2, 1=>I2) * put(1, 1=>X+im*Y)  # |0⟩⟨1|
    S1_12 = put(2, 1=>X+im*Y) * put(1, 1=>I2)  # |1⟩⟨2|
    
    # 第二个量子比特的跃迁算符
    S2_01 = put(2, 2=>I2) * put(1, 2=>X+im*Y)  # |0⟩⟨1|
    S2_12 = put(2, 2=>X+im*Y) * put(1, 2=>I2)  # |1⟩⟨2|
    
    return S1_01, S1_12, S2_01, S2_12
end

# 创建含时哈密顿量
function three_level_hamiltonian_time_dependent(Ω1_func, Ω2_func, Ω3_func, Ω4_func, Δ1, Δ2, Δ3, Δ4)
    P1_1, P1_2, P1_3, P2_1, P2_2, P2_3 = create_projectors()
    S1_01, S1_12, S2_01, S2_12 = create_transition_operators()
    
    function H(t)
        Ω1 = Ω1_func(t)  # 第一个量子比特的0-1跃迁
        Ω2 = Ω2_func(t)  # 第一个量子比特的1-2跃迁
        Ω3 = Ω3_func(t)  # 第二个量子比特的0-1跃迁
        Ω4 = Ω4_func(t)  # 第二个量子比特的1-2跃迁
        
        # 构建哈密顿量
        H = Ω1/2 * (S1_01 + S1_01') +  # 第一个量子比特的0-1耦合
            Ω2/2 * (S1_12 + S1_12') +  # 第一个量子比特的1-2耦合
            Ω3/2 * (S2_01 + S2_01') +  # 第二个量子比特的0-1耦合
            Ω4/2 * (S2_12 + S2_12') +  # 第二个量子比特的1-2耦合
            Δ1 * P1_2 +               # 第一个量子比特的|1⟩能级失谐
            Δ2 * P1_3 +               # 第一个量子比特的|2⟩能级失谐
            Δ3 * P2_2 +               # 第二个量子比特的|1⟩能级失谐
            Δ4 * P2_3                 # 第二个量子比特的|2⟩能级失谐
        return H
    end
    return H
end

# 演化函数
function evolve_three_level(reg, hamiltonian, t_total::Float64, Nt=10000)
    dt = t_total / Nt
    times = Float64[]
    populations = Vector{Float64}[]
    P1_1, P1_2, P1_3, P2_1, P2_2, P2_3 = create_projectors()
    
    # 初始布居数
    push!(times, 0.0)
    push!(populations, [
        expect(P1_1, reg),  # 第一个量子比特的|0⟩态
        expect(P1_2, reg),  # 第一个量子比特的|1⟩态
        expect(P1_3, reg),  # 第一个量子比特的|2⟩态
        expect(P2_1, reg),  # 第二个量子比特的|0⟩态
        expect(P2_2, reg),  # 第二个量子比特的|1⟩态
        expect(P2_3, reg)   # 第二个量子比特的|2⟩态
    ])
    
    # 时间演化
    for it = 1:Nt
        t = it * dt
        h = hamiltonian(t)
        apply!(reg, time_evolve(h, dt))
        
        # 记录时间和布居数
        push!(times, t)
        push!(populations, [
            expect(P1_1, reg),
            expect(P1_2, reg),
            expect(P1_3, reg),
            expect(P2_1, reg),
            expect(P2_2, reg),
            expect(P2_3, reg)
        ])
    end
    
    return reg, times, populations
end

# 创建初始态
function create_initial_state()
    # 初始态为两个量子比特都在基态 |00⟩
    reg = zero_state(2)
    return reg
end

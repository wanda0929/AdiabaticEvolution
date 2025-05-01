using Yao, Yao.EasyBuild
#using KrylovKit
#import KrylovKit: eigsolve
#using Optimisers
#using CairoMakie
using LinearAlgebra
#import LinearAlgebra: norm

# Import expect function from Yao
import Yao: expect

function Rabi_ARP_2(Ω::Float64, τ::Float64, α::Float64)
    function O(t::Float64)
        if -4 <= t <= 4
            return Ω * (exp(-(t)^2/(2*τ^2)))
        else
            return 0.0
        end
    end
    return O
end

function detuning_ARP_2(α::Float64)
    function δ(t::Float64)
        return α * t
    end
    return δ
end

function ARP_hamiltonian_2atoms_2(Ω::Float64, τ::Float64, α::Float64, U::Float64)
    h_xx = kron(X,I2) + kron(I2,X)  
    h_z = kron(-Z+I2,I2) + kron(I2,-Z+I2)
    h_i = kron(-Z+I2, -Z+I2)
    h(t) = Rabi_ARP_2(Ω, τ, α)(t)/2 * put(2,(1,2)=>h_xx) - detuning_ARP_2(α)(t)/2 * put(2,(1,2)=>h_z) + U * put(2,(1,2)=>h_i)
    #h(t) = Rabi_ARP(Ω, τ, α)(t)/2 * put(2,1=>X) + detuning_ARP(α)(t) * put(2,1=>Z) + Rabi_ARP(Ω, τ, α)(t)/2 * put(2,2=>X) + detuning_ARP(α)(t) * put(2,2=>Z) + U * put(2,(1,2)=>h_i)
    return h
end

function ARP_hamiltonian_singleatom_2(Ω::Float64, τ::Float64, α::Float64)
    h_z = -Z+I2 
    h(t) = Rabi_ARP_2(Ω, τ, α)(t)/2 * put(1,1=>X) - detuning_ARP_2(α)(t)/2 * put(1,1=>h_z) 
    return h
end

function evolve_ARP_2(reg, hamiltonian, t_total::Float64, Nt = 10000)
    dt = t_total / Nt
    phases = Float64[]
    times = Float64[]
    populations = Float64[]
    populations1 = Float64[]  # 第一个原子的布居数
    populations2 = Float64[]  # 第二个原子的布居数
    
    # 定义测量算符
    zz = (I2-Z)/2
    z2 = kron(I2,zz) + kron(zz,I2)
    
    # 初始相位和布居数
    push!(phases, 0.0)
    push!(times, 0.0)
    
    # 计算初始布居数
    if length(reg.state) == 2
        # 单量子比特情况
        push!(populations, expect(zz, reg))
    elseif length(reg.state) == 4
        # 双量子比特情况
        push!(populations, expect(z2, reg))
       #push!(populations1, expect(zz, reg |> focus!(1)))
        #push!(populations2, expect(zz, reg |> focus!(2)))
    else
        error("reg must be a single-qubit or two-qubit state")
    end
    
    # 计算一个完整周期后的演化算符矩阵
    for it = 1:Nt
        t = (it-0.5) * dt
        h = hamiltonian(t)
        apply!(reg, time_evolve(h, dt))
        
        # 计算当前状态的相位
        current_phase = angle(reg.state[1])
        push!(phases, current_phase)
        push!(times, t)
        
        # 计算激发态布居数
        if length(reg.state) == 2
            # 单量子比特情况
            excited_pop = expect(zz, reg)
            push!(populations, excited_pop)
        elseif length(reg.state) == 4
            # 双量子比特情况
            excited_pop = expect(z2, reg)
            push!(populations, excited_pop)
            #pop1 = expect(zz, reg |> focus!(1))
            #push!(populations1, pop1)
            #pop2 = expect(zz, reg |> focus!(2))
            #push!(populations2, pop2)
        else
            error("reg must be a single-qubit or two-qubit state")
        end
    end
    
    #if length(reg.state) == 2
        return reg, times, phases, populations
    #else
     #   return reg, times, phases, populations
    #end
end

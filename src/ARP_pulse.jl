using Yao, Yao.EasyBuild
#using KrylovKit
#import KrylovKit: eigsolve
#using Optimisers
#using CairoMakie
using LinearAlgebra
#import LinearAlgebra: norm

function Rabi_ARP(Ω::Float64, τ::Float64, α::Float64)
    function O(t::Float64)
        if -8 <= t <= 0
            return Ω * (exp(-(t+5)^2/(2*τ^2)))
            #return real(Ω * (exp(-(t+5)^2/(2*τ^2)-im * α * (t+5)^2/2)))
        elseif 0 < t <= 8
            return real(Ω * (exp(-(t-5)^2/(2*τ^2)))) 
            #return -real(Ω * (exp(-(t-5)^2/(2*τ^2)-im * α * (t-5)^2/2)))
        else
            return 0.0
        end
    end
    return O
end


function detuning_ARP(α::Float64, t0 = [-5.0 , 5.0])
    function δ(t::Float64)
        if -8 <= t <= 0
            return α * (t - t0[1])
        elseif 0 < t <= 8
            return α * (t - t0[2])
        else
            return 0.0
        end
    end
    return δ
end

function ARP_hamiltonian_2atoms(Ω::Float64, τ::Float64, α::Float64, U::Float64)
    h_xx = kron(X,I2) + kron(I2,X)  
    h_z = kron(-Z+I2,I2) + kron(I2,-Z+I2)
    h_i = kron(-Z+I2, -Z+I2)
    h(t) = Rabi_ARP(Ω, τ, α)(t)/2 * put(2,(1,2)=>h_xx) - detuning_ARP(α)(t)/2 * put(2,(1,2)=>h_z) + U * put(2,(1,2)=>h_i)
    #h(t) = Rabi_ARP(Ω, τ, α)(t)/2 * put(2,1=>X) + detuning_ARP(α)(t) * put(2,1=>Z) + Rabi_ARP(Ω, τ, α)(t)/2 * put(2,2=>X) + detuning_ARP(α)(t) * put(2,2=>Z) + U * put(2,(1,2)=>h_i)
    return h
end

function ARP_hamiltonian_singleatom(Ω::Float64, τ::Float64, α::Float64)
    h_z = -Z+I2 
    h(t) = Rabi_ARP(Ω, τ, α)(t)/2 * put(1,1=>X) - detuning_ARP(α)(t)/2 * put(1,1=>h_z) 
    return h
end

function evolve_ARP(reg, hamiltonian, t_total::Float64, Nt = 10000)
    dt = t_total / Nt
    phases = Float64[]
    times = Float64[]
    statevectorr = Float64[]
    zz_1 = (I2-Z)/2
    zz_2 = (I2+Z)/2
    z1 = kron(I2,zz_1) + kron(zz_1,I2)
    z2 = kron(zz_2,zz_2)
    populations_g = Float64[]
    populations_r = Float64[]
    #statevec = []
    # 初始相位
    push!(phases, 0.0)
    push!(times, 0.0)
    push!(statevectorr, 0.0)
    #push!(populations_g, 1.0)
    #push!(populations_r, 0.0)

    if length(reg.state) == 2
        push!(populations_g, expect(put(1,1=>zz_1), reg))
        push!(populations_r, expect(put(1,1=>zz_2), reg))
    elseif length(reg.state) == 4   
        push!(populations_g, expect(put(2,(1,2)=>z2), reg))
        push!(populations_r, expect(put(2,(1,2)=>z1), reg))
    else
        error("reg must be a single-qubit or two-qubit state")
    end
    #push!(populations, [abs2(reg.state[1]), abs2(reg.state[2]), abs2(reg.state[3])])
    
    # 计算一个完整周期后的演化算符矩阵
    for it = 1:Nt
        t = (it-0.5) * dt
        h = hamiltonian(t-8)
        apply!(reg, time_evolve(h, dt))
        
        # 计算当前状态的相位
        current_phase = angle(reg.state[1])
        push!(phases, current_phase)
        push!(times, t)
        if length(reg.state) == 2
            push!(populations_g, expect(put(1,1=>zz_1), reg))
            push!(populations_r, expect(put(1,1=>zz_2), reg))
            statevec1 = real(reg.state[1])
            push!(statevectorr, statevec1)
        elseif length(reg.state) == 4   
            push!(populations_g, expect(put(2,(1,2)=>z2), reg))
            push!(populations_r, expect(put(2,(1,2)=>z1), reg))
            statevec0 = real(reg.state[1])
            push!(statevectorr, statevec0)
        end
    end
    return reg, times, phases, populations_g, populations_r, statevectorr
end

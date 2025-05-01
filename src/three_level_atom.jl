using QuantumOptics
using LinearAlgebra

# Define the basis for three-level system
b = NLevelBasis(3)

# Define the operators
σ12 = transition(b, 1, 2)  # |1⟩⟨2|
σ21 = transition(b, 2, 1)  # |2⟩⟨1|
σ23 = transition(b, 2, 3)  # |2⟩⟨3|
σ32 = transition(b, 3, 2)  # |3⟩⟨2|
σ22 = transition(b, 2, 2)  # |2⟩⟨2|
σ33 = transition(b, 3, 3)  # |3⟩⟨3|

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

function three_level_hamiltonian_time_dependent(Ω1_func, Ω2_func, Δ1, Δ2)
    function H(t::Float64, ψ::Ket)
        # Create the Hamiltonian
        H_op = Ω1_func(t)/2 * (σ12 + σ21) +  # Rabi coupling between |1⟩ and |2⟩
               Ω2_func(t)/2 * (σ23 + σ32) +  # Rabi coupling between |2⟩ and |3⟩
               Δ1 * σ22 +                     # Detuning of |2⟩
               Δ2 * σ33                       # Detuning of |3⟩
        return H_op
    end
    return H
end

function evolve_three_level(ψ0, H_func, t_total, Nt=1000)
    # Create time points
    tspan = [0.0:t_total/Nt:t_total;]
    
    # Define the output function to calculate populations
    #function fout(t, ψ::Ket)
        # Calculate populations using the state vector
        #p1 = abs2(ψ.data[1])
        #p2 = abs2(ψ.data[2])
        #p3 = abs2(ψ.data[1])
        #return [p1, p2, p3]
    #end
    
    # Evolve the state using schroedinger_dynamic
    tout, ψt = timeevolution.schroedinger_dynamic(tspan, ψ0, H_func)
    #populations = [[abs2(ψt[i].data[1]), abs2(ψt[i].data[2]), abs2(ψt[i].data[3])] for i in 1:length(ψt)]
    populations = [[abs2(ψ.data[1]), abs2(ψ.data[2]), abs2(ψ.data[3])] for ψ in ψt]
    
    
    return ψt[end], tout, populations
end

# Helper function to create initial state
function create_initial_state(level::Int)
    ψ0 = basisstate(b, level)
    return ψ0
end 
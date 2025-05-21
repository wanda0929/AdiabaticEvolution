using Yao
using LinearAlgebra

n_qutrits = 3
ham_matrix = rand(ComplexF64, 3^n_qutrits, 3^n_qutrits)
ham_matrix .+= ham_matrix'

ham_3level = GeneralMatrixBlock(ham_matrix;nlevel=3)
dt = 2 * pi
time_evo_op = time_evolve(ham_3level,dt)

@assert mat(time_evo_op)  * mat(time_evo_op)'  ≈ Matrix{ComplexF64}(I, 3^n_qutrits, 3^n_qutrits)

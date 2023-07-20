# Run from the top directory of this project

using Oscar 

currentDir = pwd() 
include(joinpath(currentDir, "src/fileHandling.jl"))
include(joinpath(currentDir, "src/matroid_realization.jl"))
include(joinpath(currentDir, "src/Jacobian_Criterion.jl"))

multivariate_principal = vec(readlines("d4n9/data/multivariate_principal_ideal_4_9.dat"))

test_multivariate = []
for cir_Qstr in multivariate_principal
    A = [i for i in 1:9 if string(cir_Qstr[i]) == "*"]     
    Qstr = cir_Qstr[10:135]
    Q = matroid_from_revlex_basis_encoding(Qstr, 4, 9) 
    MR = new_matroid_realization_space(Q, A; F=QQ, saturate=true)
        
    MR = reduce_ideal_full(MR)
    I = MR.defining_ideal
    length_gens = length(gens(I))
    length_vs = length(ideal_vars(gens(I))) 
    
    x = (length(Igens) == 1 && length(ideal_vars(Igens))>1)
    
    push!(test_multivariate,x)
       
end

println(all(test_multivariate))


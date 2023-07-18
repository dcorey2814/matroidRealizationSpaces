# Run from the top directory of this project

using Oscar 

currentDir = pwd() 
include(joinpath(currentDir, "src/fileHandling.jl"))
include(joinpath(currentDir, "src/matroid_realization.jl"))
include(joinpath(currentDir, "src/Jacobian_Criterion.jl"))

univariate_ideal = vec(readlines(joinpath(currentDir, "d4n9/univariate_ideal_4_9.dat")));

test_univariate = []
for cir_Qstr in univariate_ideal
    A = [i for i in 1:9 if string(cir_Qstr[i]) == "*"]     
    Qstr = cir_Qstr[10:135]
    Q = matroid_from_revlex_basis_encoding(Qstr, 4, 9) 
    MR = new_matroid_realization_space(Q, A; F=QQ, saturate=true)
        
    MR = reduce_ideal_full(MR)
    I = MR.defining_ideal
    length_vs = length(ideal_vars(gens(I))) 
    push!(test_univariate, isone(length_vs))
end

println(all(test_univariate))

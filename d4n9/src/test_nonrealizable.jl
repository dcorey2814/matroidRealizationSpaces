# Run from the top directory of this project

using Oscar 

currentDir = pwd() 
include(joinpath(currentDir, "src/fileHandling.jl"))
include(joinpath(currentDir, "src/matroid_realization.jl"))
include(joinpath(currentDir, "src/Jacobian_Criterion.jl"))


not_realizable = vec(readlines(joinpath(currentDir, "d4n9/nonrealizable_4_9.dat")));

test_nonrealizable = []

for Qstr in not_realizable
    A = [i for i in 1:9 if string(cir_Qstr[i]) == "*"]     
    Qstr = cir_Qstr[10:135]
    Q = matroid_from_revlex_basis_encoding(Qstr, 4, 9) 
    MR = new_matroid_realization_space(Q, As[1]; F=QQ, saturate=true)
    push!(test_nonrealizable, !MR.representable)
end

println(all(test_nonrealizable))

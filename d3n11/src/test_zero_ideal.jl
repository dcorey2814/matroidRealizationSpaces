# Run from the top directory of this project

using Oscar 

currentDir = pwd() 
include(joinpath(currentDir, "src/fileHandling.jl"))
include(joinpath(currentDir, "src/isolate3Lines.jl"))
include(joinpath(currentDir, "src/matroid_realization.jl"))
include(joinpath(currentDir, "src/Jacobian_Criterion.jl"))

zero_ideal = vec(readlines(joinpath(currentDir, "d3n11/data/zero_ideal_3_11.dat")));

test_zero = []
for cir_Qstr in zero_ideal
    A = [i for i in 1:11 if string(cir_Qstr[i]) == "*"]     
    Qstr = cir_Qstr[12:176]
    Q = matroid_from_revlex_basis_encoding(Qstr, 3, 11) 
    MR = new_matroid_realization_space(Q, A; F=QQ, saturate=true)
        
    MR = reduce_ideal_full(MR)
    I = MR.defining_ideal
    push!(test_zero, iszero(I))
end


println(all(test_zero))


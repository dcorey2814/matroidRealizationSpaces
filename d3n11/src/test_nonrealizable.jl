# Run from the top directory of this project

using Oscar 

currentDir = pwd() 
include(joinpath(currentDir, "src/fileHandling.jl"))
include(joinpath(currentDir, "src/isolate3Lines.jl"))
include(joinpath(currentDir, "src/matroid_realization.jl"))
include(joinpath(currentDir, "src/Jacobian_Criterion.jl"))


not_realizable = vec(readlines(joinpath(currentDir, "d3n11/nonrealizable.dat")));

test_nonrealizable = []

for Qstr in not_realizable
    Q = matroid_from_revlex_basis_encoding(Qstr, 3, 11) 
    As = rank_plus1_circuits(Q) 
    MR = new_matroid_realization_space(Q, As[1]; F=QQ, saturate=true)
    push!(test_nonrealizable, !MR.representable)
end

println(all(test_nonrealizable))

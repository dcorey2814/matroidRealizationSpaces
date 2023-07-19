using Oscar 

currentDir = pwd() # make sure you are running from the main directory.
include(joinpath(currentDir, "src/fileHandling.jl"))
include(joinpath(currentDir, "src/isolate3Lines.jl"))
include(joinpath(currentDir, "src/matroid_realization.jl"))
include(joinpath(currentDir, "src/Jacobian_Criterion.jl"))

principal_ideal = vec(readlines(joinpath(currentDir, "d3n11/data/multivariate_principal_ideal_3_11.dat")));

test_principal = []
for cir_Qstr in principal_ideal
    A = [i for i in 1:11 if string(cir_Qstr[i]) == "*"]     
    Qstr = cir_Qstr[12:176]
    Q = matroid_from_revlex_basis_encoding(Qstr, 3, 11) 
    MR = new_matroid_realization_space(Q, A; F=QQ, saturate=true)
        
    MR = reduce_ideal_full(MR)
    R = MR.ambient_ring
    x = gens(R)
    I = MR.defining_ideal
    Igens = gens(I)
    length(Igens) != 1 && error("not principal") 
    JM = jacobian_matrix(R, x, Igens)
    nr, nc = size(JM) 
    J = I + ideal(R, [JM[1,c] for c in 1:nc])
    Sing = stepwise_saturation(J, MR.inequations)
    #Sing = realization_space_2_singular_locus(MR)
    push!(test_principal, isone(Sing))
end

print(all(test_principal))

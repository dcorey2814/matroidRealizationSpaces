# Run from the top directory of this project

using Oscar 

currentDir = pwd() 
include(joinpath(currentDir, "src/fileHandling.jl"))
include(joinpath(currentDir, "src/isolate3Lines.jl"))
include(joinpath(currentDir, "src/matroid_realization.jl"))
include(joinpath(currentDir, "src/Jacobian_Criterion.jl"))

univariate_ideal = vec(readlines(joinpath(currentDir, "d3n11/data/univariate_ideal_3_11.dat")));

test_univariate = []
for cir_Qstr in univariate_ideal
    A = [i for i in 1:11 if string(cir_Qstr[i]) == "*"]     
    Qstr = cir_Qstr[12:176]
    Q = matroid_from_revlex_basis_encoding(Qstr, 3, 11) 
    MR = new_matroid_realization_space(Q, A; F=QQ, saturate=true)
        
    MR = reduce_ideal_full(MR)
    I = MR.defining_ideal
    vs = ideal_vars(gens(I))
    length(vs) != 1 && error("not univariate")
    f = gens(I)[1]; #println("f = ", f); 
    df = derivative(f, vs[1]); #println("df = ", df); 
    push!(test_univariate, isone(ideal(MR.ambient_ring, [f,df])))
end

println(all(test_univariate))

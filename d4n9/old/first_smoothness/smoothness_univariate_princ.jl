using Oscar
using Combinatorics
pm = Polymake

currentDir = pwd()

d4n9 = vec(readlines("../data/d4n9.dat"));

include("../src/matroid_realization.jl");
include("../src/reduceIdealsv2.jl")
include("../src/JacobianCriterion.jl")
include("../src/TSC_CoordRingV2.jl")

#check principal univariate
puni = vec(readlines("sorting/principal_univariate.dat"))
io1puni = open(joinpath(pwd(), "2_princ_uni_singular_4_9.dat"), "w")
io2puni = open(joinpath(pwd(), "2_princ_uni_smooth_4_9.dat"), "w")

for t in 1:length(puni)
    Mt = matroid_from_revlex_basis_encoding(puni[t],4,9)
    I = matroid_to_reduced_expression(Mt,QQ)
   
    J = simplified_2_singular_locus(I[1], I[2])
        if !(is_one(J))
            open("2_princ_uni_singular_4_9.dat", "a") do file
                write(file, String(puni[t]),"\n")
            end
        
        else
        
            open("2_princ_uni_smooth_4_9.dat", "a") do file
                write(file, String(puni[t]),"\n")
            end
         end
end
close(io1puni)
close(io2puni)



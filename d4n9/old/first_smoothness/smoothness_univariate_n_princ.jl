using Oscar
using Combinatorics
pm = Polymake

currentDir = pwd()

d4n9 = vec(readlines("../data/d4n9.dat"));

include("../src/matroid_realization.jl");
include("../src/reduceIdealsv2.jl")
include("../src/JacobianCriterion.jl")

#check not principal univariate
npuni = vec(readlines("sorting/not_principal_univariate.dat"))  
io1npuni = open(joinpath(pwd(), "2_nprinc_uni_singular_4_9.dat"), "w")
io2npuni = open(joinpath(pwd(), "2_nprinc_uni_smooth_4_9.dat"), "w")

for t in 1:length(npuni)
    Mt = matroid_from_revlex_basis_encoding(npuni[t],4,9)
    I = matroid_to_reduced_expression(Mt,QQ)
   
    J = simplified_2_singular_locus(I[1], I[2])
        if !(is_one(J))
        
            open("2_nprinc_uni_singular_4_9.dat", "a") do file
                write(file, String(npuni[t]),"\n")
            end
        else
             open("2_nprinc_uni_smooth_4_9.dat", "a") do file
                 write(file, String(npuni[t]),"\n")
            end
         end
        
end
close(io1npuni)
close(io2npuni)


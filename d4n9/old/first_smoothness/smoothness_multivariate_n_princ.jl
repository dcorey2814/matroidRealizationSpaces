using Oscar
using Combinatorics
pm = Polymake

currentDir = pwd()

d4n9 = vec(readlines("../data/d4n9.dat"));

include("../src/matroid_realization.jl");
include("../src/reduceIdealsv2.jl")
include("../src/JacobianCriterion.jl")


#check not principal multivariate
npmult = vec(readlines("sorting/not_principal_multivariate.dat"))
io1npmult = open(joinpath(pwd(), "nprinc_mult_singular_4_9.dat"), "w")
io2npmult = open(joinpath(pwd(), "nprinc_mult_smooth_4_9.dat"), "w")

for t in 1:length(npmult)
    Mt = matroid_from_revlex_basis_encoding(npmult[t],4,9)
    I = matroid_to_reduced_expression(Mt,QQ)
   
    J = simplified_2_singular_locus(I[1], I[2])
        if !(is_one(J))
            open("nprinc_mult_singular_4_9.dat", "a") do file
                write(file, String(npmult[t]),"\n")
            end
        else
        
            open("nprinc_mult_smooth_4_9.dat", "a") do file
                 write(file, String(npmult[t]),"\n")
            end
        end
   end

close(io1npmult)
close(io2npmult)

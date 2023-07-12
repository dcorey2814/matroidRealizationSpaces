using Oscar
using Combinatorics
pm = Polymake

currentDir = pwd()

d4n9 = vec(readlines("../data/d4n9.dat"));

include("../src/matroid_realization.jl");
include("../src/reduceIdealsv2.jl")
include("../src/JacobianCriterion.jl")
include("../src/TSC_CoordRingV2.jl")

#check principal multivariate

pmult = vec(readlines("sorting/principal_multivariate.dat"))
io1pmult = open(joinpath(pwd(), "2_princ_mult_singular_4_9.dat"), "w")
io2pmult = open(joinpath(pwd(), "2_princ_mult_smooth_4_9.dat"), "w")

for t in 1:length(pmult)
    Mt = matroid_from_revlex_basis_encoding(pmult[t],4,9)
    I = matroid_to_reduced_expression(Mt,QQ)
   
    J = simplified_2_singular_locus(I[1], I[2])
        if !(is_one(J))
            open("2_princ_mult_singular_4_9.dat", "a") do file
                write(file, String(pmult[t]),"\n")
            end
        else
            open("2_princ_mult_smooth_4_9.dat", "a") do file
                 write(file, String(pmult[t]),"\n")
            end
        end
        
end
close(io1pmult)
close(io2pmult)


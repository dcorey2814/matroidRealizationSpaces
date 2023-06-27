using Oscar
using Combinatorics
pm = Polymake

currentDir = pwd()

d4n8 = vec(readlines("../data/d4n8.dat"));

include("../src/matroid_realization.jl");
include("../src/reduceIdealsv2.jl")
include("../src/JacobianCriterion.jl")
include("../src/TSC_CoordRingV2.jl")

weird = vec(readlines("sorting/weird.dat"))
sp_48 = load("4planes_4_8.dat")


matroids_4_8 = union(
#vec(readlines("d3n11/sorting/int_generator.dat")),
vec(readlines("sorting/principal_multivariate.dat")),
vec(readlines("sorting/principal_univariate.dat")),
vec(readlines("sorting/not_principal_multivariate.dat")),
vec(readlines("sorting/not_principal_univariate.dat"))    
)

io1 = open(joinpath(pwd(), "singular_4_8_2.dat"), "w")
io2 = open(joinpath(pwd(), "smooth_4_8_2.dat"), "w")

for t in 1:length(matroids_4_8)
    Mt = matroid_from_revlex_basis_encoding(matroids_4_8[t],4,8)
    I = matroid_to_reduced_expression(Mt,QQ)
   
    J = simplified_2_singular_locus(I[1], I[2])
        if !(is_one(J))
            write(io1, String(matroids_4_8[t]),"\n")
           #close(io1)
            #continue
        else
             write(io2, String(matroids_4_8[t]),"\n")
             #close(io2)
             #continue
         end
        
end
close(io1)
close(io2)


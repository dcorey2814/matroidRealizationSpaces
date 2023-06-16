using Oscar
using Combinatorics
pm = Polymake

include("../src/matroid_realization.jl");
include("../src/reduceIdealsv2.jl")
include("../src/JacobianCriterion.jl")

matroids_3_9 = union(
#vec(readlines("d3n9/sorting/int_generator.dat")),
vec(readlines("d3n9/sorting/principal_multivariate.dat")),
vec(readlines("d3n9/sorting/principal_univariate.dat")),
#vec(readlines("d3n9/sorting/not_principal_multivariate.dat"))
)

io1 = open(joinpath(pwd(), "d3n9/check_singular_3_9.dat"), "w")
io2 = open(joinpath(pwd(), "d3n9/catch_3_9.dat"), "w")

for t in 1:length(matroids_3_9)
    Mt = matroid_from_revlex_basis_encoding(matroids_3_9[t],3,9)
    I = matroid_to_reduced_expression(Mt,QQ)
   
    J = simplified_2_singular_locus(I[1], I[2])
        if !(is_one(J))
            write(io1, String(matroids_3_9[t]),"\n")
           #close(io1)
            #continue
        else
             write(io2, String(matroids_3_9[t]),"\n")
             #close(io2)
             #continue
         end
        
end
close(io1)
close(io2)

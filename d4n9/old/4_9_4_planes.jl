using Oscar
using Combinatorics
pm = Polymake


include("../src/matroid_realization.jl");
include("../src/reduceIdealsv2.jl")
include("../src/JacobianCriterion.jl")

d4n9 = vec(readlines("../data/d4n9.dat"))

println("check function")
function in_planes(n,L)

    Ln = [l for l in L if n in l]

    return(length(Ln))

end

println("isolate matroids with such that every element of [1,9] is in at least 4 planes")

for t in 1:length(d4n9)
    
    println(t)

    Mt = matroid_from_revlex_basis_encoding(d4n9[t],4,9)

    L = [h for h in hyperplanes(Mt) if length(h)>3]

    ns = [n for n in 1:9 if in_planes(n,L)>3]

    if length(ns) == 9
                println("four planes")
                    open("4_9_four_planes.dat", "a") do file
                         write(file,d4n9[t],"\n")
                    end 
    end
end
    

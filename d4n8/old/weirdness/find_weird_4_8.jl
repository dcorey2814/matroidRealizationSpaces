using Oscar
using Combinatorics
pm = Polymake


include("../src/matroid_realization.jl");
include("../src/reduceIdealsv2.jl")
include("../src/JacobianCriterion.jl")

d4n8 = vec(readlines("../data/d4n8.dat"))

println("check function")


println("find 4-8 matroids with no reference circuit")

for t in 1:length(d4n8)
    
    println(t)

    Mt = matroid_from_revlex_basis_encoding(d4n8[t],4,8)

    C = [c for c in circuits(Mt) if length(c) == rank(Mt)+1]

  
    if length(C) == 0
                println("wierd")
                    open("4_8_weird.dat", "a") do file
                         write(file,d4n8[t],"\n")
                    end 
    end
end

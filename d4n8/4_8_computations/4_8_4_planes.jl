using Oscar
using Combinatorics
pm = Polymake



d4n8 = vec(readlines("d4n8.dat"))

println("check function")
function in_planes(n,L)

    Ln = [l for l in L if n in l]

    return(length(Ln))

end

println("isolate matroids with such that every element of [1,8] is in at least 4 planes")

for t in 1:length(d4n8)
    
    println(t)

    Mt = matroid_from_revlex_basis_encoding(d4n8[t],4,8)

    L = [h for h in hyperplanes(Mt) if length(h)>3]

    ns = [n for n in 1:8 if in_planes(n,L)>3]


    if (length(ns) == 8 && is_connected(Mt))
                println("four planes")
                    open("4_8_four_planes_connected.dat", "a") do file
                         write(file,d4n8[t],"\n")
                    end 
    end
end
    
